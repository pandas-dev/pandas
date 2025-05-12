window.addEventListener("DOMContentLoaded", function() {
  var absBaseUrl = document.baseURI;
  var baseUrl = location.protocol + "//" + location.hostname
  if (location.port) {
     baseUrl = baseUrl + ":" + location.port
  }
  var currentLanguage = document.documentElement.lang;
  var languages = JSON.parse(document.getElementById("languages").getAttribute('data-lang').replace(/'/g, '"'));
  const languageNames = {
        'en': 'English',
        'es': 'Español',
        'fr': 'Français',
        'pt': 'Português'
      }

  // Handle preview URLs on github
  // If preview URL changes, this regex will need to be updated
  const re = /preview\/pandas-dev\/pandas\/(?<pr>[0-9]*)\//g;
  var previewUrl = '';
  for (const match of absBaseUrl.matchAll(re)) {
    previewUrl = `/preview/pandas-dev/pandas/${match.groups.pr}`;
  }
  var pathName = location.pathname.replace(previewUrl, '')

  // Create dropdown menu
  function makeDropdown(options) {
    var dropdown = document.createElement("li");
    dropdown.classList.add("nav-item");
    dropdown.classList.add("dropdown");

    var link = document.createElement("a");
    link.classList.add("nav-link");
    link.classList.add("dropdown-toggle");
    link.setAttribute("data-bs-toggle", "dropdown");
    link.setAttribute("href", "#");
    link.setAttribute("role", "button");
    link.setAttribute("aria-haspopup", "true");
    link.setAttribute("aria-expanded", "false");
    link.textContent = languageNames[currentLanguage];

    var dropdownMenu = document.createElement("div");
    dropdownMenu.classList.add("dropdown-menu");

    options.forEach(function(i) {
      var dropdownItem = document.createElement("a");
      dropdownItem.classList.add("dropdown-item");
      dropdownItem.textContent = languageNames[i] || i.toUpperCase();
      dropdownItem.setAttribute("href", "#");
      dropdownItem.addEventListener("click", function() {
        var urlLanguage = '';
        if (i !== 'en') {
          urlLanguage = '/' + i;
        }
        pathName = pathName.replace('/' + currentLanguage + '/', '/')
        var newUrl = baseUrl + previewUrl + urlLanguage + pathName
        window.location.href = newUrl;
      });
      dropdownMenu.appendChild(dropdownItem);
    });

    dropdown.appendChild(link);
    dropdown.appendChild(dropdownMenu);
    return dropdown;
  }

  var container = document.getElementById("language-switcher-container");
  if (container) {
    var dropdown = makeDropdown(languages);
    container.appendChild(dropdown);
  }
});
