import { BlobFile } from "generic-filehandle";
import { vizsnp } from "./vizsnp";

const resultsElement = document.querySelector("#results");
const vizsnpFormElement = document.querySelector("#vizsnpForm");

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
    resultsElement.innerHTML = "Please upload a VCF and TBI file";
    return;
  }

  try {
    // show loading message and disable submit button
    resultsElement.innerHTML = "Loading...";
    form.submit.disabled = true;

    // run vizsnp
    let start = Date.now();
    const vcfFileInput = new BlobFile(vcfFile);
    const tbiFileInput = new BlobFile(tbiFile);
    const results = await vizsnp(vcfFileInput, tbiFileInput, limit, gene, "human");
    let end = Date.now();

    // show results
    resultsElement.innerHTML = `Loaded ${results.length} variants in ${(end - start) / 1000}s\n`;
    resultsElement.innerHTML += JSON.stringify(results, null, 2);
  } catch (e) {
    console.error(e);
    resultsElement.innerHTML = e.toString();
  } finally {
    // re-enable submit button
    form.submit.disabled = false;
  }
};
