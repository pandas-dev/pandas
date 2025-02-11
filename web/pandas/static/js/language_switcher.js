window.addEventListener("DOMContentLoaded", function() {

  // ABS_BASE_URL is the full URL
  // e.g. https://pandas.pydata.org/en/getting_started.html
  var ABS_BASE_URL = document.baseURI;
  var CURRENT_LANGUAGE = document.documentElement.lang;
  var languages = JSON.parse(document.getElementById("languages").getAttribute('data-lang').replace(/'/g, '"'));

  const language_names = {
        'en': 'English',
        'es': 'Español',
        'zh': '中文',
        'ja': '日本語',
        'ko': '한국어',
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
    link.textContent = language_names[CURRENT_LANGUAGE];

    var dropdownMenu = document.createElement("div");
    dropdownMenu.classList.add("dropdown-menu");

    options.forEach(function(i) {
      var dropdownItem = document.createElement("a");
      dropdownItem.classList.add("dropdown-item");
      dropdownItem.textContent = language_names[i] || i.toUpperCase();
      dropdownItem.setAttribute("href", "#");
      dropdownItem.addEventListener("click", function() {
        var newUrl = ABS_BASE_URL.replace('/' + CURRENT_LANGUAGE + '/', '/' + i + '/');
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
