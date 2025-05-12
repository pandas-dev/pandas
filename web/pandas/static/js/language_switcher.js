window.addEventListener("DOMContentLoaded", function() {
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
        var pathName = location.pathname.replace('/' + currentLanguage + '/', '/')
        var newUrl = baseUrl + urlLanguage + pathName
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
