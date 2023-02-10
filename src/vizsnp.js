import { TabixIndexedFile } from "@gmod/tabix";
import VCF from "@gmod/vcf";

/**
 * Returns a promise that resolves after a given number of milliseconds
 * @param {number} ms 
 * @returns 
 */
export function sleep(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

/**
 * Encodes an object as a x-www-form-urlencoded string
 * @param {any} object 
 * @returns 
 */
export function encodeUrl(object) {
  return Object.entries(object)
    .map(([k, v]) => `${k}=${Array.isArray(v) ? v.join(",") : v}`)
    .join("&");
}

/**
 * Splits an array into chunks of a given size
 * @param {any[]} array 
 * @param {number} chunkSize 
 * @returns 
 */
export function chunk(array, chunkSize = 200) {
  let chunks = [];
  for (let i = 0; i < array.length; i += chunkSize) {
    chunks.push(array.slice(i, i + chunkSize));
  }
  return chunks;
}

/**
 * Retrieves protein ids from the uniprot idmapping service
 * @param {string|string[]} ids - List of ids to map
 * @param {string} from - default: Ensembl. Map from this id type
 * @param {string} to - default: UniProtKB-Swiss-Prot. Map to this id type
 * @returns ID mapping results
 */
export async function getProteinIds(ids, from = "Ensembl", to = "UniProtKB-Swiss-Prot") {
  const url = "https://rest.uniprot.org/idmapping";

  // submit uniprot idmapping job
  const runResponse = await fetch(`${url}/run`, {
    method: "post",
    body: encodeUrl({ from, to, ids }),
    headers: {
      "Content-Type": "application/x-www-form-urlencoded",
    },
  });
  const { jobId } = await runResponse.json();

  // wait for job to finish
  let status = { results: null };
  while (!status.results) {
    const statusResponse = await fetch(`${url}/status/${jobId}`);
    status = await statusResponse.json();
    if (!status.results) {
      await sleep(1000);
    }
  }

  return status.results;
}

/**
 * Retrieves a list of PDB structures which map to each UniProt accession id, sorted by coverage of the protein and, if the same, resolution.
 * @param {string[]} uniprotAccessionIds - UniProt accession ids to get best structures for
 * @returns
 */
export async function getBestStructures(uniprotAccessionIds) {
  const url = "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures";
  const response = await fetch(url, {
    method: "post",
    body: uniprotAccessionIds.join(","),
    headers: {
      "Content-Type": "application/x-www-form-urlencoded",
    },
  });
  return await response.json();
}

/**
 * Returns variant consequences for multiple regions
 * @param {string[]} variants - List of variants to get consequences for
 * @param {string} species - default: human. Species to use for variant consequence prediction
 * @returns
 */
export async function getVariantConsequences(variants, species = "human") {
  const url = `https://rest.ensembl.org/vep/${species}/region?uniprot=1`;
  const response = await fetch(url, {
    method: "post",
    headers: {
      "Content-Type": "application/json",
    },
    body: JSON.stringify({ variants }),
  });
  return await response.json();
}

/**
 * Processes variant consequences and appends PDB structure and ICN3D links
 * @param {data} data - Variant consequence data
 * @returns 
 */
export async function processVariantConsequences(data) {
  let results = [];
  for (let result of data) {
    if (result.transcript_consequences) {
      const consequence = result.transcript_consequences.find((c) => (c.sift_prediction || c.polyphen_prediction) && c.swissprot?.length);
      if (consequence) {
        const swissprotId = consequence.swissprot[0].split(".")[0];
        const aminoAcids = consequence.amino_acids.split("/");
        const aminoAcidMutation = [aminoAcids[0], consequence.protein_start, aminoAcids[1]].join("");
        results.push({
          variant: result.input,
          geneId: consequence.gene_id,
          swissprotId,
          pdbId: null,
          aminoAcidMutation,
          siftPrediction: consequence.sift_prediction,
          siftScore: consequence.sift_score,
          polyphenPrediction: consequence.polyphen_prediction,
          polyphenScore: consequence.polyphen_score,
          icn3d: 'TODO: Generate ICN3D link',
          consequence,
        });
      }
    }
  }

  // get best structures for each protein
  const bestStructures = await getBestStructures(results.map(r => r.swissprotId));
  for (let key in bestStructures) {
    results = results.map(r => {
      if (r.swissprotId === key) {
        const pdbStructures = bestStructures[key];
        const highestResolutionStructure = pdbStructures[0];
        r.pdbId = [highestResolutionStructure.pdb_id, highestResolutionStructure.chain_id].join('_');
        r.pdbStructures = pdbStructures;
      }
      return r;
    });
  }

  // todo: get icn3d link for each variant

  return results;
}

/**
 * Reads a VCF file using tabix and returns a list of variants
 * @param {string} filepath 
 * @returns 
 */
async function readVcfFile(vcfFile, tbiFile, limit = 1000, filter = null) {
  const tbiIndexed = new TabixIndexedFile({ 
    filehandle: vcfFile,
    tbiFilehandle: tbiFile
  });
  const headerText = await tbiIndexed.getHeader();
  const tbiVCFParser = new VCF({ header: headerText });
  const chromosomes = Object.keys(tbiVCFParser.metadata.contig);

  const variants = [];
  for (let chromosome of chromosomes) {
    await tbiIndexed.getLines(chromosome, 0, undefined, (line) => {
      if (variants.length < limit && (!filter || line.includes(filter))) {
        variants.push(tbiVCFParser.parseLine(line));
      }
    });
  }
  return variants;
}

/**
 * Formats a variant for use with Ensembl VEP
 * @param {any} variant - VCF variant
 * @returns 
 */
function parseVariant(variant) {
  return [variant.CHROM.replace("chr", ""), variant.POS, variant.ID ? variant.ID[0] : ".", variant.REF, variant.ALT.join(","), "."].join(" ");
}

/**
 * Executes the VizSNP pipeline
 * @param {*} vcfFilePath 
 * @param {*} species 
 * @returns 
 */
export async function vizsnp(vcfFile, tbiFile, limit = 1000, filter = null, species = "human") {
  const vcfFileContents = await readVcfFile(vcfFile, tbiFile, limit, filter);
  const variantChunks = chunk(vcfFileContents, 200);

  let results = [];
  for (const chunk of variantChunks) {
    const variants = chunk.map(parseVariant);
    const variantConsequences = await getVariantConsequences(variants, species);
    const processedVariantConsequences = await processVariantConsequences(variantConsequences)
    results = results.concat(processedVariantConsequences);
  }

  return results;
}
