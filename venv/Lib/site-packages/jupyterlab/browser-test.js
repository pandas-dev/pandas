/*
 * Copyright (c) Jupyter Development Team.
 * Distributed under the terms of the Modified BSD License.
 */

const playwright = require('playwright');
const path = require('path');
const fs = require('fs');

// eslint-disable-next-line no-redeclare
const URL = process.argv[2];
const BROWSER_VAR = 'JLAB_BROWSER_TYPE';
const BROWSER = process.env[BROWSER_VAR] || 'chromium';
const HEADLESS_VAR = 'JLAB_BROWSER_HEADLESS';
const HEADLESS = process.env[HEADLESS_VAR] === 'false' ? false : true;
const OUTPUT_VAR = 'JLAB_BROWSER_CHECK_OUTPUT';
const OUTPUT = process.env[OUTPUT_VAR];

let nextScreenshot = 0;
const screenshotStem = `screenshot-${+new Date()}`;

if (OUTPUT) {
  console.log(`Screenshots will be saved in ${OUTPUT}...`);
  if (!fs.existsSync(OUTPUT)) {
    console.log(`Creating ${OUTPUT}...`);
    fs.mkdirSync(OUTPUT, { recursive: true });
  }
}

async function main() {
  /* eslint-disable no-console */
  console.info(`Starting headless ${BROWSER}...`);
  let testError = null;

  const pwBrowser = playwright[BROWSER];
  const browser = await pwBrowser.launch({
    headless: HEADLESS,
    logger: {
      isEnabled: () => !!OUTPUT,
      log: (name, severity, message, args) => console.log(name, message)
    }
  });

  const context = await browser.newContext();
  const page = await context.newPage();

  async function screenshot() {
    if (!OUTPUT) {
      return;
    }
    const screenshotPath = path.join(
      OUTPUT,
      `${screenshotStem}-${++nextScreenshot}.png`
    );
    console.log(`Capturing screenshot ${screenshotPath}...`);
    await page.screenshot({
      type: 'png',
      path: screenshotPath
    });
  }

  console.info('Navigating to page:', URL);
  await page.goto(URL);
  console.info('Waiting for page to load...');

  // Wait for the local file to redirect on notebook >= 6.0
  await page.waitForNavigation();

  console.log('Waiting for page content..');

  try {
    await page.locator('#jupyter-config-data').waitFor({ state: 'attached' });
  } catch (reason) {
    console.error('Error loading JupyterLab page:', reason);
    // Limit to 1000 characters
    console.error((await page.content()).substring(0, 1000));
    testError = reason;
  }

  console.log('Waiting for #main selector...');
  await page.waitForSelector('#main', { timeout: 100000 });

  console.log('Waiting for #browserTest selector...');
  const el = await page.waitForSelector('#browserTest', {
    timeout: 100000,
    state: 'attached'
  });
  console.log('Waiting for application to start...');

  try {
    await page.waitForSelector('.completed', { state: 'attached' });
  } catch (e) {
    testError = e;
  }

  await screenshot();

  const textContent = await el.getProperty('textContent');
  const errors = JSON.parse(await textContent.jsonValue());

  for (let error of errors) {
    console.error(`Parsed an error from text content: ${error.message}`, error);
    testError = true;
  }

  await browser.close();

  if (testError) {
    throw testError;
  }
  console.info('Browser test complete');
}

// Stop the process if an error is raised in the async function.
process.on('unhandledRejection', up => {
  throw up;
});

void main();
