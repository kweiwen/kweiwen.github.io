const STORAGE_KEY = "lang";
const CONTENT_VERSION = "1";
const REVEAL_TIMEOUT_MS = 1500;

const documentRoot = document.documentElement;
const bioList = document.querySelector("[data-i18n-list]");
const bioItemTemplate = bioList?.firstElementChild?.cloneNode(true);

const dictionaries = new Map();

let activeLanguage = documentRoot.dataset.lang || "en";
let activeDictionary;
let onLanguageApplied = () => {};
let revealTimeoutId = window.setTimeout(revealContent, REVEAL_TIMEOUT_MS);

function revealContent() {
  window.clearTimeout(revealTimeoutId);
  delete documentRoot.dataset.i18nPending;
}

function resolve(dictionary, path) {
  return path.split(".").reduce((value, key) => (value == null ? value : value[key]), dictionary);
}

async function loadDictionary(language) {
  if (!dictionaries.has(language)) {
    const response = await fetch(`content/${language}.json?v=${CONTENT_VERSION}`);

    if (!response.ok) {
      throw new Error(`Missing content for "${language}"`);
    }

    dictionaries.set(language, await response.json());
  }

  return dictionaries.get(language);
}

function renderBioList(items) {
  if (!bioList || !bioItemTemplate || !Array.isArray(items)) {
    return;
  }

  const fragment = document.createDocumentFragment();

  for (const { date, text } of items) {
    const item = bioItemTemplate.cloneNode(true);
    item.querySelector("time").textContent = date;
    item.querySelector("p").textContent = text;
    fragment.append(item);
  }

  bioList.replaceChildren(fragment);
}

function applyDictionary(dictionary) {
  activeDictionary = dictionary;

  document.querySelectorAll("[data-i18n]").forEach((element) => {
    const value = resolve(dictionary, element.dataset.i18n);

    if (typeof value === "string") {
      element.textContent = value;
    }
  });

  document.querySelectorAll("[data-i18n-attr]").forEach((element) => {
    for (const binding of element.dataset.i18nAttr.split(";")) {
      const [attribute, key] = binding.split(":");
      const value = resolve(dictionary, key.trim());

      if (typeof value === "string") {
        element.setAttribute(attribute.trim(), value);
      }
    }
  });

  renderBioList(resolve(dictionary, bioList?.dataset.i18nList ?? ""));

  documentRoot.lang = dictionary.htmlLang ?? activeLanguage;
  documentRoot.dataset.lang = activeLanguage;
  revealContent();
  onLanguageApplied();
}

export function t(key) {
  const value = resolve(activeDictionary, key);
  return typeof value === "string" ? value : undefined;
}

export function getLanguage() {
  return activeLanguage;
}

export async function setLanguage(language) {
  if (language === activeLanguage && activeDictionary) {
    return;
  }

  const previousLanguage = activeLanguage;
  activeLanguage = language;

  try {
    const dictionary = await loadDictionary(language);
    localStorage.setItem(STORAGE_KEY, language);
    applyDictionary(dictionary);
  } catch {
    activeLanguage = previousLanguage;
    revealContent();
  }
}

export function toggleLanguage() {
  return setLanguage(activeLanguage === "zh" ? "en" : "zh");
}

export function initI18n(onApply) {
  onLanguageApplied = onApply ?? onLanguageApplied;
  return setLanguage(activeLanguage);
}
