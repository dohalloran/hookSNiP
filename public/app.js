
(function(){
  "use strict";
  const $  = (s) => document.querySelector(s);
  const $$ = (s) => document.querySelectorAll(s);

  function wireUI(){
    const syncGenome = () => {
      const mode = document.querySelector("input[name='srcmode']:checked").value;
      $("#genomeSelect").disabled = (mode !== "built-in");
      $("#customGenome").disabled = (mode !== "custom");
    };
    $$("input[name='srcmode']").forEach(r => r.addEventListener("change", syncGenome));
    syncGenome();
    $("#runBtn").addEventListener("click", runDesign);
  }

  async function runDesign(){
    $("#results").innerHTML = "";
    $("#downloadBox").classList.add("hidden");
    const status = $("#runStatus"); const prog = $("#runProgress");
    prog.hidden = false; prog.value = 0; status.textContent = "Uploading files...";

    try{
      const mode = document.querySelector("input[name='srcmode']:checked").value;
      const flank = Math.max(50, parseInt($("#flank").value,10) || 300);
      const top = Math.max(1, parseInt($("#topPerAllele").value,10) || 1);
      const mism = $("#insertMismatch").checked;

      const fd = new FormData();
      fd.append("srcmode", mode);
      fd.append("flank", String(flank));
      fd.append("top_per_allele", String(top));
      fd.append("insert_mismatch", String(mism));

      if (mode === "built-in") {
        const name = $("#genomeSelect").value;
        fd.append("genome_name", name);
      } else {
        const gf = $("#customGenome").files[0];
        if (!gf) throw new Error("Please choose a custom FASTA (.fa/.fasta/.fna or .fa.gz).");
        fd.append("genome", gf, gf.name);
      }
      const gff = $("#gffFile").files[0];
      if (!gff) throw new Error("Please upload a GFF/GFF3 file.");
      fd.append("gff3", gff, gff.name);

      const vcf = $("#vcfFile").files[0];
      if (!vcf) throw new Error("Please upload a VCF file.");
      fd.append("vcf", vcf, vcf.name);

      status.textContent = "Server is analyzing files...";
      const resp = await fetch("/api/run", { method: "POST", body: fd });
      if (!resp.ok) {
        const t = await resp.text();
        throw new Error("Server error: " + resp.status + " " + t);
      }
      const out = await resp.json();
      if (!out.ok) throw new Error(out.error || "Unknown error");

      $("#results").innerHTML = out.table_html || "(no table)";
      const blob = new Blob([out.primers_html || ""], { type: "text/html" });
      const url = URL.createObjectURL(blob);
      $("#downloadHtmlBtn").onclick = () => {
        const a=document.createElement("a");
        a.href = url; a.download = "primers.html"; a.click();
      };
      $("#downloadBox").classList.remove("hidden");
      status.textContent = "Done.";
    } catch (e){
      console.error(e);
      status.textContent = "Error: " + (e && e.message ? e.message : e);
    } finally {
      prog.hidden = true;
    }
  }

  document.addEventListener("DOMContentLoaded", wireUI);
})();
