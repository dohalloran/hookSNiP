(function(){
  "use strict";
  const $  = (s) => document.querySelector(s);
  const $$ = (s) => document.querySelectorAll(s);

  function wireUI(){
    const syncGenome = () => {
      const mode = document.querySelector("input[name='srcmode']:checked").value;
      $("#genomeSelect").disabled  = (mode !== "built-in");
      $("#customGenome").disabled  = (mode !== "custom");
    };
    $$("input[name='srcmode']").forEach(r => r.addEventListener("change", syncGenome));
    syncGenome();

    // FIX #11: The original wired the click handler as addEventListener("click", runDesign).
    // Because runDesign is async, any unhandled rejection (e.g. a network error before
    // the try-block's catch runs) would silently disappear. Wrap in an explicit .catch().
    $("#runBtn").addEventListener("click", () => {
      runDesign().catch(e => {
        const status = $("#runStatus");
        status.textContent = "Unexpected error: " + (e && e.message ? e.message : String(e));
        status.classList.add("error");
        $("#runProgress").hidden = true;
        console.error(e);
      });
    });
  }

  async function runDesign(){
    $("#results").innerHTML = "";
    $("#downloadBox").classList.add("hidden");
    const status = $("#runStatus");
    const prog   = $("#runProgress");

    // FIX #12: Clear error class from any previous failed run so success
    // renders correctly and doesn't stay red.
    status.classList.remove("error");

    prog.hidden = false;
    prog.value  = 0;
    status.textContent = "Uploading files...";

    try {
      const mode  = document.querySelector("input[name='srcmode']:checked").value;
      const flank = Math.max(50, parseInt($("#flank").value, 10) || 300);
      const top   = Math.max(1,  parseInt($("#topPerAllele").value, 10) || 1);
      const mism  = $("#insertMismatch").checked;

      const fd = new FormData();
      fd.append("srcmode",        mode);
      fd.append("flank",          String(flank));
      fd.append("top_per_allele", String(top));
      fd.append("insert_mismatch", String(mism));

      if (mode === "built-in") {
        fd.append("genome_name", $("#genomeSelect").value);
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

      // FIX #13: Add an AbortController so requests don't hang forever if the
      // local server goes away mid-run. 10-minute timeout is generous for large genomes.
      const controller = new AbortController();
      const timeoutId  = setTimeout(() => controller.abort(), 10 * 60 * 1000);

      let resp;
      try {
        resp = await fetch("/api/run", { method: "POST", body: fd, signal: controller.signal });
      } finally {
        clearTimeout(timeoutId);
      }

      if (!resp.ok) {
        const t = await resp.text().catch(() => "");
        let detail = "";
        try { detail = JSON.parse(t).detail || t; } catch { detail = t; }
        throw new Error(`Server error ${resp.status}: ${detail}`);
      }

      const out = await resp.json();
      if (!out.ok) throw new Error(out.error || "Unknown server error");

      $("#results").innerHTML = out.table_html || "(no table)";

      // FIX #14: Revoke the previous object URL before creating a new one
      // to avoid accumulating blob URLs across multiple runs.
      if (window._primersBlobUrl) {
        URL.revokeObjectURL(window._primersBlobUrl);
      }
      const blob = new Blob([out.primers_html || ""], { type: "text/html" });
      window._primersBlobUrl = URL.createObjectURL(blob);
      $("#downloadHtmlBtn").onclick = () => {
        const a = document.createElement("a");
        a.href = window._primersBlobUrl;
        a.download = "primers.html";
        a.click();
      };

      $("#downloadBox").classList.remove("hidden");
      const processed = out.info && out.info.processed != null ? out.info.processed : "?";
      status.textContent = `Done. Variants processed: ${processed}.`;

    } catch (e) {
      console.error(e);
      // FIX #15: Distinguish abort (timeout) from other errors for a clearer message.
      const msg = e.name === "AbortError"
        ? "Request timed out (server took > 10 min). Try a smaller input or increase timeout."
        : (e && e.message ? e.message : String(e));
      status.textContent = "Error: " + msg;
      status.classList.add("error");
    } finally {
      prog.hidden = true;
    }
  }

  document.addEventListener("DOMContentLoaded", wireUI);
})();
