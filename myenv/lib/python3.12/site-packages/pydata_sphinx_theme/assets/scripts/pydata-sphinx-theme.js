// Define the custom behavior of the page
import { documentReady } from "./mixin";
import { compare, validate } from "compare-versions";

import "../styles/pydata-sphinx-theme.scss";

/*******************************************************************************
 * Theme interaction
 */

var prefersDark = window.matchMedia("(prefers-color-scheme: dark)");

/**
 * set the the body theme to the one specified by the user browser
 *
 * @param {event} e
 */
function autoTheme(e) {
  document.documentElement.dataset.theme = prefersDark.matches
    ? "dark"
    : "light";
}

/**
 * Set the theme using the specified mode.
 * It can be one of ["auto", "dark", "light"]
 *
 * @param {str} mode
 */
function setTheme(mode) {
  if (mode !== "light" && mode !== "dark" && mode !== "auto") {
    console.error(`Got invalid theme mode: ${mode}. Resetting to auto.`);
    mode = "auto";
  }

  // get the theme
  var colorScheme = prefersDark.matches ? "dark" : "light";
  document.documentElement.dataset.mode = mode;
  var theme = mode == "auto" ? colorScheme : mode;
  document.documentElement.dataset.theme = theme;
  // TODO: remove this line after Bootstrap upgrade
  // v5.3 has a colors mode: https://getbootstrap.com/docs/5.3/customize/color-modes/
  document.querySelectorAll(".dropdown-menu").forEach((el) => {
    if (theme === "dark") {
      el.classList.add("dropdown-menu-dark");
    } else {
      el.classList.remove("dropdown-menu-dark");
    }
  });

  // save mode and theme
  localStorage.setItem("mode", mode);
  localStorage.setItem("theme", theme);
  console.log(`[PST]: Changed to ${mode} mode using the ${theme} theme.`);

  // add a listener if set on auto
  prefersDark.onchange = mode == "auto" ? autoTheme : "";
}

/**
 * Change the theme option order so that clicking on the btn is always a change
 * from "auto"
 */
function cycleMode() {
  const defaultMode = document.documentElement.dataset.defaultMode || "auto";
  const currentMode = localStorage.getItem("mode") || defaultMode;

  var loopArray = (arr, current) => {
    var nextPosition = arr.indexOf(current) + 1;
    if (nextPosition === arr.length) {
      nextPosition = 0;
    }
    return arr[nextPosition];
  };

  // make sure the next theme after auto is always a change
  var modeList = prefersDark.matches
    ? ["auto", "light", "dark"]
    : ["auto", "dark", "light"];
  var newMode = loopArray(modeList, currentMode);
  setTheme(newMode);
}

/**
 * add the theme listener on the btns of the navbar
 */
function addModeListener() {
  // the theme was set a first time using the initial mini-script
  // running setMode will ensure the use of the dark mode if auto is selected
  setTheme(document.documentElement.dataset.mode);

  // Attach event handlers for toggling themes colors
  document.querySelectorAll(".theme-switch-button").forEach((el) => {
    el.addEventListener("click", cycleMode);
  });
}

/*******************************************************************************
 * TOC interactivity
 */

/**
 * TOC sidebar - add "active" class to parent list
 *
 * Bootstrap's scrollspy adds the active class to the <a> link,
 * but for the automatic collapsing we need this on the parent list item.
 *
 * The event is triggered on "window" (and not the nav item as documented),
 * see https://github.com/twbs/bootstrap/issues/20086
 */
function addTOCInteractivity() {
  window.addEventListener("activate.bs.scrollspy", function () {
    const navLinks = document.querySelectorAll(".bd-toc-nav a");

    navLinks.forEach((navLink) => {
      navLink.parentElement.classList.remove("active");
    });

    const activeNavLinks = document.querySelectorAll(".bd-toc-nav a.active");
    activeNavLinks.forEach((navLink) => {
      navLink.parentElement.classList.add("active");
    });
  });
}

/*******************************************************************************
 * Scroll
 */

/**
 * Navigation sidebar scrolling to active page
 */
function scrollToActive() {
  // If the docs nav doesn't exist, do nothing (e.g., on search page)
  if (!document.querySelector(".bd-docs-nav")) {
    return;
  }

  var sidebar = document.querySelector("div.bd-sidebar");

  // Remember the sidebar scroll position between page loads
  // Inspired on source of revealjs.com
  let storedScrollTop = parseInt(
    sessionStorage.getItem("sidebar-scroll-top"),
    10
  );

  if (!isNaN(storedScrollTop)) {
    // If we've got a saved scroll position, just use that
    sidebar.scrollTop = storedScrollTop;
    console.log("[PST]: Scrolled sidebar using stored browser position...");
  } else {
    // Otherwise, calculate a position to scroll to based on the lowest `active` link
    var sidebarNav = document.querySelector(".bd-docs-nav");
    var active_pages = sidebarNav.querySelectorAll(".active");
    if (active_pages.length > 0) {
      // Use the last active page as the offset since it's the page we're on
      var latest_active = active_pages[active_pages.length - 1];
      var offset =
        latest_active.getBoundingClientRect().y -
        sidebar.getBoundingClientRect().y;
      // Only scroll the navbar if the active link is lower than 50% of the page
      if (latest_active.getBoundingClientRect().y > window.innerHeight * 0.5) {
        let buffer = 0.25; // Buffer so we have some space above the scrolled item
        sidebar.scrollTop = offset - sidebar.clientHeight * buffer;
        console.log("[PST]: Scrolled sidebar using last active link...");
      }
    }
  }

  // Store the sidebar scroll position
  window.addEventListener("beforeunload", () => {
    sessionStorage.setItem("sidebar-scroll-top", sidebar.scrollTop);
  });
}

/*******************************************************************************
 * Search
 */

/**
 * Find any search forms on the page and return their input element
 */
var findSearchInput = () => {
  let forms = document.querySelectorAll("form.bd-search");
  if (!forms.length) {
    // no search form found
    return;
  } else {
    var form;
    if (forms.length == 1) {
      // there is exactly one search form (persistent or hidden)
      form = forms[0];
    } else {
      // must be at least one persistent form, use the first persistent one
      form = document.querySelector(
        "div:not(.search-button__search-container) > form.bd-search"
      );
    }
    return form.querySelector("input");
  }
};

/**
 * Activate the search field on the page.
 * - If there is a search field already visible it will be activated.
 * - If not, then a search field will pop up.
 */
var toggleSearchField = () => {
  // Find the search input to highlight
  let input = findSearchInput();

  // if the input field is the hidden one (the one associated with the
  // search button) then toggle the button state (to show/hide the field)
  let searchPopupWrapper = document.querySelector(".search-button__wrapper");
  let hiddenInput = searchPopupWrapper.querySelector("input");
  if (input === hiddenInput) {
    searchPopupWrapper.classList.toggle("show");
  }
  // when toggling off the search field, remove its focus
  if (document.activeElement === input) {
    input.blur();
  } else {
    input.focus();
    input.select();
    input.scrollIntoView({ block: "center" });
  }
};

/**
 * Add an event listener for toggleSearchField() for Ctrl/Cmd + K
 */
var addEventListenerForSearchKeyboard = () => {
  window.addEventListener(
    "keydown",
    (event) => {
      let input = findSearchInput();
      // toggle on Ctrl+k or ⌘+k
      if ((event.ctrlKey || event.metaKey) && event.code == "KeyK") {
        event.preventDefault();
        toggleSearchField();
      }
      // also allow Escape key to hide (but not show) the dynamic search field
      else if (document.activeElement === input && event.code == "Escape") {
        toggleSearchField();
      }
    },
    true
  );
};

/**
 * Change the search hint to `meta key` if we are a Mac
 */
var changeSearchShortcutKey = () => {
  let forms = document.querySelectorAll("form.bd-search");
  var isMac = window.navigator.platform.toUpperCase().indexOf("MAC") >= 0;
  if (isMac) {
    forms.forEach(
      (f) => (f.querySelector("kbd.kbd-shortcut__modifier").innerText = "⌘")
    );
  }
};

/**
 * Activate callbacks for search button popup
 */
var setupSearchButtons = () => {
  changeSearchShortcutKey();
  addEventListenerForSearchKeyboard();

  // Add the search button trigger event callback
  document.querySelectorAll(".search-button__button").forEach((btn) => {
    btn.onclick = toggleSearchField;
  });

  // Add the search button overlay event callback
  let overlay = document.querySelector(".search-button__overlay");
  if (overlay) {
    overlay.onclick = toggleSearchField;
  }
};

/*******************************************************************************
 * Version Switcher
 * Note that this depends on two variables existing that are defined in
 * and `html-page-context` hook:
 *
 * - DOCUMENTATION_OPTIONS.pagename
 * - DOCUMENTATION_OPTIONS.theme_switcher_url
 */

/**
 * Check if corresponding page path exists in other version of docs
 * and, if so, go there instead of the homepage of the other docs version
 *
 * @param {event} event the event that trigger the check
 */
async function checkPageExistsAndRedirect(event) {
  // ensure we don't follow the initial link
  event.preventDefault();
  let currentFilePath = `${DOCUMENTATION_OPTIONS.pagename}.html`;
  let tryUrl = event.currentTarget.getAttribute("href");
  let otherDocsHomepage = tryUrl.replace(currentFilePath, "");
  try {
    let head = await fetch(tryUrl, { method: "HEAD" });
    if (head.ok) {
      location.href = tryUrl; // the page exists, go there
    } else {
      location.href = otherDocsHomepage;
    }
  } catch (err) {
    // something went wrong, probably CORS restriction, fallback to other docs homepage
    location.href = otherDocsHomepage;
  }
}

/**
 * Load and parse the version switcher JSON file from an absolute or relative URL.
 *
 * @param {string} url The URL to load version switcher entries from.
 */
async function fetchVersionSwitcherJSON(url) {
  // first check if it's a valid URL
  try {
    var result = new URL(url);
  } catch (err) {
    if (err instanceof TypeError) {
      // assume we got a relative path, and fix accordingly. But first, we need to
      // use `fetch()` to follow redirects so we get the correct final base URL
      const origin = await fetch(window.location.origin, { method: "HEAD" });
      result = new URL(url, origin.url);
    } else {
      // something unexpected happened
      throw err;
    }
  }
  // load and return the JSON
  const response = await fetch(result);
  const data = await response.json();
  return data;
}

// Populate the version switcher from the JSON data
function populateVersionSwitcher(data, versionSwitcherBtns) {
  const currentFilePath = `${DOCUMENTATION_OPTIONS.pagename}.html`;
  versionSwitcherBtns.forEach((btn) => {
    // Set empty strings by default so that these attributes exist and can be used in CSS selectors
    btn.dataset["activeVersionName"] = "";
    btn.dataset["activeVersion"] = "";
  });
  // in case there are multiple entries with the same version string, this helps us
  // decide which entry's `name` to put on the button itself. Without this, it would
  // always be the *last* version-matching entry; now it will be either the
  // version-matching entry that is also marked as `"preferred": true`, or if that
  // doesn't exist: the *first* version-matching entry.
  data = data.map((entry) => {
    // does this entry match the version that we're currently building/viewing?
    entry.match =
      entry.version == DOCUMENTATION_OPTIONS.theme_switcher_version_match;
    entry.preferred = entry.preferred || false;
    // if no custom name specified (e.g., "latest"), use version string
    if (!("name" in entry)) {
      entry.name = entry.version;
    }
    return entry;
  });
  const hasMatchingPreferredEntry = data
    .map((entry) => entry.preferred && entry.match)
    .some(Boolean);
  var foundMatch = false;
  // create links to the corresponding page in the other docs versions
  data.forEach((entry) => {
    // create the node
    const anchor = document.createElement("a");
    anchor.setAttribute(
      "class",
      "dropdown-item list-group-item list-group-item-action py-1"
    );
    anchor.setAttribute("href", `${entry.url}${currentFilePath}`);
    anchor.setAttribute("role", "option");
    const span = document.createElement("span");
    span.textContent = `${entry.name}`;
    anchor.appendChild(span);
    // Add dataset values for the version and name in case people want
    // to apply CSS styling based on this information.
    anchor.dataset["versionName"] = entry.name;
    anchor.dataset["version"] = entry.version;
    // replace dropdown button text with the preferred display name of the
    // currently-viewed version, rather than using sphinx's {{ version }} variable.
    // also highlight the dropdown entry for the currently-viewed version's entry
    let matchesAndIsPreferred = hasMatchingPreferredEntry && entry.preferred;
    let matchesAndIsFirst =
      !hasMatchingPreferredEntry && !foundMatch && entry.match;
    if (matchesAndIsPreferred || matchesAndIsFirst) {
      anchor.classList.add("active");
      versionSwitcherBtns.forEach((btn) => {
        btn.innerText = entry.name;
        btn.dataset["activeVersionName"] = entry.name;
        btn.dataset["activeVersion"] = entry.version;
      });
      foundMatch = true;
    }
    // There may be multiple version-switcher elements, e.g. one
    // in a slide-over panel displayed on smaller screens.
    document.querySelectorAll(".version-switcher__menu").forEach((menu) => {
      // we need to clone the node for each menu, but onclick attributes are not
      // preserved by `.cloneNode()` so we add onclick here after cloning.
      let node = anchor.cloneNode(true);
      node.onclick = checkPageExistsAndRedirect;
      // on click, AJAX calls will check if the linked page exists before
      // trying to redirect, and if not, will redirect to the homepage
      // for that version of the docs.
      menu.append(node);
    });
  });
}

/*******************************************************************************
 * Warning banner when viewing non-stable version of the docs.
 */

/**
 * Show a warning banner when viewing a non-stable version of the docs.
 *
 * adapted 2023-06 from https://mne.tools/versionwarning.js, which was
 * originally adapted 2020-05 from https://scikit-learn.org/versionwarning.js
 *
 * @param {Array} data The version data used to populate the switcher menu.
 */
function showVersionWarningBanner(data) {
  const version = DOCUMENTATION_OPTIONS.theme_version;
  // figure out what latest stable version is
  var preferredEntries = data.filter((entry) => entry.preferred);
  if (preferredEntries.length !== 1) {
    const howMany = preferredEntries.length == 0 ? "No" : "Multiple";
    console.log(
      `[PST] ${howMany} versions marked "preferred" found in versions JSON, ignoring.`
    );
    return;
  }
  const preferredVersion = preferredEntries[0].version;
  const preferredURL = preferredEntries[0].url;
  // if already on preferred version, nothing to do
  const versionsAreComparable = validate(version) && validate(preferredVersion);
  if (versionsAreComparable && compare(version, preferredVersion, "=")) {
    return;
  }
  // now construct the warning banner
  var outer = document.createElement("div");
  const middle = document.createElement("div");
  const inner = document.createElement("div");
  const bold = document.createElement("strong");
  const button = document.createElement("a");
  // these classes exist since pydata-sphinx-theme v0.10.0
  outer.classList = "bd-header-version-warning container-fluid";
  middle.classList = "bd-header-announcement__content";
  inner.classList = "sidebar-message";
  button.classList =
    "sd-btn sd-btn-danger sd-shadow-sm sd-text-wrap font-weight-bold ms-3 my-1 align-baseline";
  button.href = `${preferredURL}${DOCUMENTATION_OPTIONS.pagename}.html`;
  button.innerText = "Switch to stable version";
  button.onclick = checkPageExistsAndRedirect;
  // add the version-dependent text
  inner.innerText = "This is documentation for ";
  const isDev =
    version.includes("dev") ||
    version.includes("rc") ||
    version.includes("pre");
  const newerThanPreferred =
    versionsAreComparable && compare(version, preferredVersion, ">");
  if (isDev || newerThanPreferred) {
    bold.innerText = "an unstable development version";
  } else if (versionsAreComparable && compare(version, preferredVersion, "<")) {
    bold.innerText = `an old version (${version})`;
  } else {
    bold.innerText = `version ${version}`;
  }
  outer.appendChild(middle);
  middle.appendChild(inner);
  inner.appendChild(bold);
  inner.appendChild(document.createTextNode("."));
  inner.appendChild(button);
  document.body.prepend(outer);
}

/*******************************************************************************
 * MutationObserver to move the ReadTheDocs button
 */

/**
 * intercept the RTD flyout and place it in the rtd-footer-container if existing
 * if not it stays where on top of the page
 */
function initRTDObserver() {
  const mutatedCallback = (mutationList, observer) => {
    mutationList.forEach((mutation) => {
      // Check whether the mutation is for RTD, which will have a specific structure
      if (mutation.addedNodes.length === 0) {
        return;
      }
      if (mutation.addedNodes[0].data === undefined) {
        return;
      }
      if (mutation.addedNodes[0].data.search("Inserted RTD Footer") != -1) {
        mutation.addedNodes.forEach((node) => {
          document.getElementById("rtd-footer-container").append(node);
        });
      }
    });
  };

  const observer = new MutationObserver(mutatedCallback);
  const config = { childList: true };
  observer.observe(document.body, config);
}

// fetch the JSON version data (only once), then use it to populate the version
// switcher and maybe show the version warning bar
var versionSwitcherBtns = document.querySelectorAll(
  ".version-switcher__button"
);
const hasSwitcherMenu = versionSwitcherBtns.length > 0;
const hasVersionsJSON = DOCUMENTATION_OPTIONS.hasOwnProperty(
  "theme_switcher_json_url"
);
const wantsWarningBanner = DOCUMENTATION_OPTIONS.show_version_warning_banner;

if (hasVersionsJSON && (hasSwitcherMenu || wantsWarningBanner)) {
  const data = await fetchVersionSwitcherJSON(
    DOCUMENTATION_OPTIONS.theme_switcher_json_url
  );
  populateVersionSwitcher(data, versionSwitcherBtns);
  if (wantsWarningBanner) {
    showVersionWarningBanner(data);
  }
}

/*******************************************************************************
 * Call functions after document loading.
 */

documentReady(addModeListener);
documentReady(scrollToActive);
documentReady(addTOCInteractivity);
documentReady(setupSearchButtons);
documentReady(initRTDObserver);
