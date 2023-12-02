var sd_labels_by_text = {};

function ready() {
  const li = document.getElementsByClassName("sd-tab-label");
  for (const label of li) {
    syncId = label.getAttribute("data-sync-id");
    if (syncId) {
      label.onclick = onLabelClick;
      if (!sd_labels_by_text[syncId]) {
        sd_labels_by_text[syncId] = [];
      }
      sd_labels_by_text[syncId].push(label);
    }
  }
}

function onLabelClick() {
  // Activate other inputs with the same sync id.
  syncId = this.getAttribute("data-sync-id");
  for (label of sd_labels_by_text[syncId]) {
    if (label === this) continue;
    label.previousElementSibling.checked = true;
  }
  window.localStorage.setItem("sphinx-design-last-tab", syncId);
}

document.addEventListener("DOMContentLoaded", ready, false);
