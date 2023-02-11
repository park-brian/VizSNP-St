import { BlobFile } from "generic-filehandle";
import { vizsnp } from "./vizsnp";

const resultsElement = document.querySelector("#results");
const resultsGridElement = document.querySelector("#resultsGrid");
const vizsnpFormElement = document.querySelector("#vizsnpForm");
const loadingElement = document.querySelector("#loading");

vizsnpFormElement.onsubmit = async function submit(event) {
  event.preventDefault();

  // get form values
  const form = event.target;
  const inputFiles = form.inputFiles.files;
  const limit = +form.limit.value || 1000;
  const gene = form.gene.value;

  // get vcf and tbi files
  const vcfFile = [...inputFiles].find((file) => file.name.endsWith(".vcf.gz") || file.name.endsWith(".vcf"));
  const tbiFile = [...inputFiles].find((file) => file.name.endsWith(".tbi"));

  // validate input
  if (!vcfFile || !tbiFile) {
    resultsElement.innerText = "Please upload a VCF and TBI file.";
    return;
  }

  try {
    // show loading message and disable submit button
    resultsElement.innerHTML = "";
    resultsGridElement.innerHTML = "";
    form.submit.disabled = true;
    loadingElement.hidden = false;

    // run vizsnp
    let start = Date.now();
    const vcfFileInput = new BlobFile(vcfFile);
    const tbiFileInput = new BlobFile(tbiFile);
    const results = await vizsnp(vcfFileInput, tbiFileInput, limit, gene, "human");
    let end = Date.now();

    // show results in console
    resultsElement.innerText = `Loaded ${results.length} variant(s) in ${(end - start) / 1000}s. View browser console to inspect results.`;
    console.table(results);

    // link cell renderer
    function renderLink({ value }) {
      const a = document.createElement("a");
      a.href = value;
      a.target = "_blank";
      a.innerText = value;
      return a;
    }

    // show results in table
    new agGrid.Grid(resultsGridElement, {
      rowData: results,
      defaultColDef: {
        resizable: true,
      },
      columnDefs: [
        { field: "variant", headerName: "Variant" },
        { field: "geneId", headerName: "Gene" },
        { field: "uniprotId", headerName: "UniProt ID" },
        { field: "pdbId", headerName: "PDB ID" },
        { field: "aminoAcidMutation", headerName: "Amino Acid Mutation" },
        { field: "siftPrediction", headerName: "SIFT Prediction" },
        { field: "siftScore", headerName: "SIFT Score" },
        { field: "polyphenPrediction", headerName: "PolyPhen Prediction" },
        { field: "polyphenScore", headerName: "PolyPhen Score" },
        { field: "icn3dAlphafold", headerName: "iCn3D AlphaFold URL", cellRenderer: renderLink },
        { field: "icn3dPdb", headerName: "iCn3D PDB URL", cellRenderer: renderLink },
      ],
      onFirstDataRendered: ({ columnApi }) => setTimeout(() => columnApi.autoSizeAllColumns(), 0),
    });
  } catch (e) {
    console.error(e);
    resultsElement.innerText = e.toString();
  } finally {
    // re-enable submit button
    form.submit.disabled = false;
    loadingElement.hidden = true;
  }
};
