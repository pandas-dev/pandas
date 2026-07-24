/**
 * Moto Dashboard JavaScript Library
 * Handles dynamic table generation for Moto AWS service mocking data
 */

// Wait for DOM to be fully loaded before initializing the dashboard
window.addEventListener('DOMContentLoaded', function () {
  refreshData();
});

/**
 * Fetches data from the Moto API and updates the dashboard tables
 * @async
 * @function refreshData
 */
const refreshData = async () => {
  await fetch('/moto-api/data.json')
    .then(result => result.json())
    .then(data => {
      // Update the dashboard with fetched data
      updateTable(data);

      // Clean up CSS classes after initial load - hide non-active tabs
      $('.onbootonly').each(function () {
        this.className = 'tab-pane';
      });
    });
};

/**
 * Expected data structure from /moto-api/data.json:
 * {
 *   // Service name (e.g., "cloudwatch", "core", "events")
 *   "cloudwatch": {
 *     // Resource type name with array of resource instances
 *     "Alarm": [],           // Array of alarm objects
 *     "Dashboard": [],       // Array of dashboard objects
 *     "InsightRule": [],     // Array of insight rule objects
 *     "MetricAggregatedDatum": [],
 *     "MetricDatum": [],
 *     "MetricDatumBase": []
 *   },
 *   "core": {
 *     "CloudFormationModel": []
 *   },
 *   "events": {
 *     "Archive": [],
 *     "Connection": [],
 *     "Destination": [],
 *     "EventBus": [],
 *     "PartnerEventSource": [],
 *     "Replay": [],
 *     "Rule": []
 *   }
 * }
 */

/**
 * Main function to update the dashboard with service data
 * Creates tabbed interface with Bootstrap tables for each service and resource type
 * @param {Object} data - The service data object from Moto API
 * @function updateTable
 */
function updateTable(data) {
  // Get the main data container element from the HTML
  const dataContainer = document.getElementById('data');
  if (!dataContainer) {
    console.error('Element with id="data" not found');
    return;
  }

  // Create navigation menu tabs for each service
  const listService = document.getElementById('list-service');
  if (listService) {
    Object.keys(data).forEach((serviceName, svcIndex) => {
      // Create navigation item for each service
      const listItem = document.createElement('li');
      listItem.className = 'nav-item';

      // Create clickable tab button for service
      const listItemLink = document.createElement('button');
      listItemLink.className = 'nav-link';
      listItemLink.setAttribute('id', `service-tab-${serviceName}`);
      listItemLink.setAttribute('data-bs-toggle', 'pill');
      listItemLink.setAttribute('data-bs-target', `#div-table-${serviceName}`);
      listItemLink.setAttribute('role', 'tab');
      listItemLink.setAttribute('aria-controls', `div-table-${serviceName}`);
      listItemLink.setAttribute('type', 'button');

      // Set first service as active/selected by default
      if (svcIndex == 0) {
        listItemLink.setAttribute('aria-selected', 'true');
        listItemLink.classList.add('active');
        focusedService = serviceName;
      } else {
        listItemLink.setAttribute('aria-selected', 'false');
      }

      // Set service name as button text
      listItemLink.textContent = serviceName;
      listItem.appendChild(listItemLink);
      listService.appendChild(listItem);
    });
  }

  // Clear existing content from previous renders
  dataContainer.innerHTML = '';

  // Handle empty or invalid data
  if (!data || Object.keys(data).length == 0) {
    const noDataMessage = document.createElement('p');
    noDataMessage.textContent = 'No data';
    dataContainer.appendChild(noDataMessage);
    return;
  }

  // Create tab content for each service
  // Structure: Service Tab -> Resource Tables within each tab
  Object.keys(data).forEach((serviceName) => {
    const resources = data[serviceName];

    // Create container div for this service's content
    const serviceDiv = document.createElement('div');
    serviceDiv.setAttribute('id', `div-table-${serviceName}`);

    // Set CSS classes for Bootstrap tab functionality
    if (serviceName == focusedService) {
      // First service is active and visible
      serviceDiv.className = 'tab-pane show active firsttab';
    } else {
      // Other services are initially hidden (will be shown on first load, then hidden)
      serviceDiv.className = 'tab-pane show active onbootonly';
    }

    // Set ARIA attributes for accessibility
    serviceDiv.setAttribute('role', 'tabpanel');
    serviceDiv.setAttribute('aria-labelledby', `service-tab-${serviceName}`);
    serviceDiv.setAttribute('tabindex', "0");

    // Add visual separator at top of each service section
    const hr = document.createElement('hr');
    hr.className = 'border border-primary border-1 opacity-75';
    serviceDiv.appendChild(hr);

    // Create container for all resource tables within this service
    const resourcesDiv = document.createElement('div');
    resourcesDiv.setAttribute('id', `div-table-content-${serviceName}`);

    // Append resource container to service div, then service div to main container
    serviceDiv.appendChild(resourcesDiv);
    dataContainer.appendChild(serviceDiv);

    // Process each resource type within the current service
    Object.keys(resources).forEach((resourceName, index) => {
      const resourceData = resources[resourceName];
      const columns = [];

      // Add visual separator between multiple resource tables
      if (index > 0) {
        const hr = document.createElement('hr');
        resourcesDiv.appendChild(hr);
      }

      // Create title header for each resource type
      const resourceTitle = document.createElement('h4');
      resourceTitle.className = 'float-start mt-3';
      resourceTitle.textContent = `${serviceName} - ${resourceName}`;
      resourcesDiv.appendChild(resourceTitle);

      // Handle empty resource arrays - show message instead of empty table
      if (resourceData.length === 0) {
        const emptyMessage = document.createElement('p');
        emptyMessage.textContent = `No data available for ${resourceName}`;
        emptyMessage.className = 'fst-italic';

        // Remove float class since no table will be displayed
        resourceTitle.className = '';
        resourcesDiv.appendChild(emptyMessage);
        return;
      }


      // Dynamically generate table columns based on object properties
      // Collect all unique keys from all objects to handle varying object structures
      if (columns.length === 0) {
        const allKeys = new Set();

        // Scan all resource objects to find all possible property keys
        resourceData.forEach(item => {
          Object.keys(item).forEach(key => allKeys.add(key));
        });

        // Create column definition for each unique key
        allKeys.forEach(key => {
          columns.push({
            field: key,
            title: key.charAt(0).toUpperCase() + key.slice(1), // Capitalize first letter
            sortable: true
          });
        });
      }

      // Create HTML table element for Bootstrap Table
      const resourceTable = document.createElement('table');
      resourceTable.setAttribute('id', `table-${serviceName}-${resourceName}`);
      resourceTable.setAttribute('tabindex', "0");
      resourceTable.setAttribute('data-height', '500');             // Set fixed height
      resourceTable.setAttribute('data-show-columns', 'true');      // Show column visibility toggle
      resourcesDiv.appendChild(resourceTable);

      // Preprocess data: Convert complex objects to JSON strings for display
      resourceData.forEach(item => {
        Object.keys(item).forEach(key => {
          if (typeof item[key] === 'object') {
            item[key] = JSON.stringify(item[key]);
          }
        });
      });

      // Initialize Bootstrap Table with data and configuration
      $(`#table-${serviceName}-${resourceName}`).bootstrapTable({
        data: resourceData,        // The actual data array
        columns: columns,          // Column definitions
        search: true,             // Enable search functionality
      });

      // Apply Bootstrap styling classes to the table
      resourceTable.className = 'table table-striped table-bordered';
    });
  });
}


// Color mode 
/*!
 * Color mode toggler for Bootstrap's docs (https://getbootstrap.com/)
 * Copyright 2011-2025 The Bootstrap Authors
 * Licensed under the Creative Commons Attribution 3.0 Unported License.
 */

(() => {
  'use strict'

  const getStoredTheme = () => localStorage.getItem('theme')
  const setStoredTheme = theme => localStorage.setItem('theme', theme)

  const getPreferredTheme = () => {
    const storedTheme = getStoredTheme()
    if (storedTheme) {
      return storedTheme
    }

    return window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light'
  }

  const setTheme = theme => {
    if (theme === 'auto') {
      document.documentElement.setAttribute('data-bs-theme', (window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light'))
    } else {
      document.documentElement.setAttribute('data-bs-theme', theme)
    }
  }

  setTheme(getPreferredTheme())

  const showActiveTheme = (theme, focus = false) => {
    const themeSwitcher = document.querySelector('#bd-theme')

    if (!themeSwitcher) {
      return
    }

    const themeSwitcherText = document.querySelector('#bd-theme-text')
    // const activeThemeIcon = document.querySelector('.theme-icon-active use')
    const btnToActive = document.querySelector(`[data-bs-theme-value="${theme}"]`)
    // const svgOfActiveBtn = btnToActive.querySelector('svg use').getAttribute('href')

    document.querySelectorAll('[data-bs-theme-value]').forEach(element => {
      element.classList.remove('active')
      element.setAttribute('aria-pressed', 'false')
    })

    btnToActive.classList.add('active')
    btnToActive.setAttribute('aria-pressed', 'true')
    // activeThemeIcon.setAttribute('href', svgOfActiveBtn)
    const themeSwitcherLabel = `${themeSwitcherText.textContent} (${btnToActive.dataset.bsThemeValue})`
    themeSwitcher.setAttribute('aria-label', themeSwitcherLabel)

    if (focus) {
      themeSwitcher.focus()
    }
  }

  window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => {
    const storedTheme = getStoredTheme()
    if (storedTheme !== 'light' && storedTheme !== 'dark') {
      setTheme(getPreferredTheme())
    }
  })

  window.addEventListener('DOMContentLoaded', () => {
    showActiveTheme(getPreferredTheme())

    document.querySelectorAll('[data-bs-theme-value]')
      .forEach(toggle => {
        toggle.addEventListener('click', () => {
          const theme = toggle.getAttribute('data-bs-theme-value')
          setStoredTheme(theme)
          setTheme(theme)
          showActiveTheme(theme, true)
        })
      })
  })
})()