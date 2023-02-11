import { GenericFilehandle } from "generic-filehandle";
import { TabixIndexedFile } from "@gmod/tabix";
import VCF from "@gmod/vcf";

/**
 * Runs the VizSNP pipeline to retrieve variant consequences and icn3d links for alphafold/pdb structures
 * @param {GenericFilehandle} vcfFileHandle - The file handle for the VCF file
 * @param {GenericFilehandle} tbiFileHandle - The file handle for the tabix-indexed file
 * @param {number} limit - The maximum number of variants to retrieve
 * @param {string} filter - A string to filter variants by
 * @param {string} species - The species for which to retrieve variant consequences
 * @returns {Promise<any[]>} An array of variant consequences and icn3d links for alphafold/pdb structures
 */
export async function vizsnp(vcfFileHandle, tbiFileHandle, limit = 1000, filter = null, species = "human") {
  const vcfFileContents = await readVcfFile(vcfFileHandle, tbiFileHandle, limit, filter);

  let results = [];
  for await (const variants of chunk(vcfFileContents, 200)) {
    const parsedVariants = variants.map(parseVariant);
    const variantConsequences = await getVariantConsequences(parsedVariants, species);
    const processedVariantConsequences = await processVariantConsequences(variantConsequences);
    results = results.concat(processedVariantConsequences);
  }

  return results;
}

/**
 * Retrieves protein ids from the uniprot idmapping service
 * @param {string|string[]} ids - List of ids to map
 * @param {string} from - default: Ensembl. Map from this id type
 * @param {string} to - default: UniProtKB-Swiss-Prot. Map to this id type
 * @returns {any} ID mapping results
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
 * @param {string[]} uniprotIds - UniProt accession ids to get best structures for
 * @returns {Promise<any>} Best structures for each UniProt accession id
 */
export async function getBestStructures(uniprotIds) {
  const url = "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures";
  const response = await fetch(url, {
    method: "post",
    body: uniprotIds.join(","),
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
 * @returns {Promise<any[]>} Variant consequences
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
 * @returns {Promise<any[]>} Processed variant consequence data
 */
async function processVariantConsequences(data) {
  const results = [];

  for (const result of data) {
    const consequence = result?.transcript_consequences?.find((c) => (c.sift_prediction || c.polyphen_prediction) && c.swissprot?.length);
    if (!consequence) continue;

    const uniprotId = consequence.swissprot[0].split(".")[0];
    const aminoAcids = consequence.amino_acids.split("/");
    const aminoAcidMutation = [aminoAcids[0], consequence.protein_start, aminoAcids[1]].join("");

    results.push({
      variant: result.input,
      geneId: consequence.gene_id,
      uniprotId,
      pdbId: null,
      aminoAcids,
      aminoAcidMutation,
      aminoAcidMutationPosition: consequence.protein_start,
      siftPrediction: consequence.sift_prediction,
      siftScore: consequence.sift_score,
      polyphenPrediction: consequence.polyphen_prediction,
      polyphenScore: consequence.polyphen_score,
      icn3dPdb: null,
      icn3dAlphafold: null,
      consequence,
    });
  }

  // get best structures for each protein (sorted by coverage and resolution)
  const bestStructures = await getBestStructures(results.map((r) => r.uniprotId));
  for (let result of results) {
    if (!bestStructures[result.uniprotId]) continue;
    const pdbStructures = bestStructures[result.uniprotId];
    const highestResolutionStructure = pdbStructures[0];
    result.pdbId = [highestResolutionStructure.pdb_id, highestResolutionStructure.chain_id].join("_");
    result.pdbStructures = pdbStructures;
  }

  for (let result of results) {
    if (result.siftPrediction !== "deleterious" && result.polyphenPrediction !== "probably_damaging") continue;
    if (result.pdbId) result.icn3dPdb = await getIcn3dPdbUrl(result);
    result.icn3dAlphafold = await getIcn3dAlphafoldUrl(result);
  }

  return results;
}

/**
 * Retrieves ICN3D URL for an AlphaFold structure
 * @param {*} result 
 * @returns {string} ICN3D URL
 */
async function getIcn3dAlphafoldUrl(result) {
  let url = 'https://www.ncbi.nlm.nih.gov/Structure/icn3d/full.html';
  let date = new Date();
  const params = {
    afid: result.uniprotId,
    date: [date.getFullYear(), date.getMonth() + 1, date.getDate()].join(""),
    v: "3.12.7",
    command: [
      "view annotations",
      "set annotation cdd",
      "set view detailed view",
      "set thickness | stickrad 0.2",
      `add track | chainid ${result.uniprotId}_A | title SIFT_predict | text ${result.aminoAcidMutationPosition} ${result.aminoAcids[1]}`,
      `add track | chainid ${result.uniprotId}_A | title Polyphen_predict | text ${result.aminoAcidMutationPosition} ${result.aminoAcids[1]}`,
      `scap interaction ${result.uniprotId}_A_${result.aminoAcidMutationPosition}_${result.aminoAcids[1]}`
    ].join("; ")
  }
  return url + '?' + encodeUrl(params, true);
}

/**
 * Retrieves ICN3D URL for a PDB structure
 * @param {*} result 
 * @returns {string} ICN3D URL
 */
async function getIcn3dPdbUrl(result) {
  let url = 'https://www.ncbi.nlm.nih.gov/Structure/icn3d/full.html';
  let date = new Date();
  let offset = 0; // todo: get offset from pdb structure
  const params = {
    pdbid: result.pdbId,
    date: [date.getFullYear(), date.getMonth() + 1, date.getDate()].join(""),
    v: "3.12.7",
    command: [
      "view annotations",
      "set annotation cdd",
      "set view detailed view",
      "set thickness | stickrad 0.2",
      `add track | chainid ${result.pdbId} | title SIFT_predict | text ${result.aminoAcidMutationPosition + offset} ${result.aminoAcids[1]}`,
      `add track | chainid ${result.pdbId} | title Polyphen_predict | text ${result.aminoAcidMutationPosition + offset} ${result.aminoAcids[1]}`,
      `scap interaction ${result.pdbId}_${result.aminoAcidMutationPosition + offset}_${result.aminoAcids[1]}`
    ].join("; ")
  }

  return url + '?' + encodeUrl(params, true);
}

/**
 * Reads a VCF file using tabix and returns a list of variants
 * @param {string} filehandle - The filehandle to the VCF file
 * @param {string} tbiFilehandle - The filehandle to the tabix-indexed file
 * @param {number} limit - The maximum number of variants to retrieve
 * @param {string} filter - A string to filter variants by
 * @returns {Promise<any[]>} An array of variants
 */
export async function readVcfFile(filehandle, tbiFilehandle, limit = 1000, filter = null) {
  const tbiIndexed = new TabixIndexedFile({ filehandle, tbiFilehandle });
  const headerText = await tbiIndexed.getHeader();
  const tbiVCFParser = new VCF({ header: headerText });
  const refNames = Object.keys(tbiVCFParser.metadata.contig);

  const variants = [];
  for (let refName of refNames) {
    await tbiIndexed.getLines(refName, 0, undefined, (line) => {
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
  let { CHROM, POS, ID, REF, ALT } = variant;
  if (CHROM) CHROM = CHROM.replace("chr", "");
  if (ID) ID = ID[0];
  if (ALT) ALT = ALT.join(",");

  return [CHROM, POS, ID, REF, ALT, null].map((e) => e ?? ".").join(" ");
}


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
export function encodeUrl(object, encode = false) {
  return Object.entries(object)
    .map(([k, v]) => `${encode ? encodeURIComponent(k) : k}=${(Array.isArray(v) ? v : [v]).map(v => encode ? encodeURIComponent(v) : v).join(",")}`)
    .join("&");
}

/**
 * Splits an array into chunks of a given size
 * @param {any[]} array
 * @param {number} chunkSize
 * @returns
 */
export function* chunk(array, chunkSize = 200) {
  for (let i = 0; i < array.length; i += chunkSize) {
    yield array.slice(i, i + chunkSize);
  }
}
