document.addEventListener("DOMContentLoaded", function () {
  document.querySelectorAll(".copybtn").forEach((btn) => {
    btn.addEventListener("click", async function (event) {
      event.preventDefault();

      const codeBlock = btn.closest("div.highlight");
      if (!codeBlock) return;

      const clone = codeBlock.cloneNode(true);

      // Remove spans with data-hidden="true"
      clone.querySelectorAll("span.copybutton").forEach((el) => {
        if (el.getAttribute("data-hidden") === "true") {
          el.remove();
        }
      });

      const text = clone.innerText;

      try {
        await navigator.clipboard.writeText(text);
        console.log("Copied cleaned code to clipboard.");
      } catch (err) {
        console.warn("Clipboard API failed, using fallback.");
        fallbackCopy(text);
      }
    });
  });

  function fallbackCopy(text) {
    const textarea = document.createElement("textarea");
    textarea.value = text;
    document.body.appendChild(textarea);
    textarea.select();
    try {
      document.execCommand("copy");
      console.log("Copied using fallback.");
    } catch (err) {
      console.error("Fallback copy failed:", err);
    }
    document.body.removeChild(textarea);
  }
});
