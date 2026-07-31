const documentRoot = document.documentElement;
const glassRoot = document.querySelector("#glass-root");
const siteHeader = document.querySelector(".site-header");
const mobileMenuGlass = document.querySelector(".mobile-menu-glass");
const mainContent = document.querySelector("main");
const themeToggle = document.querySelector(".theme-toggle");
const menuToggle = document.querySelector(".menu-toggle");
const navigationMenu = document.querySelector(".nav-links");
const themeColorMeta = document.querySelector('meta[name="theme-color"]');
const identityContent = document.querySelector(".identity-copy");
const profileAvatar = document.querySelector("[data-profile-avatar]");

let liquidGlassInstance;
let isHeaderScrolled = false;
let scrollAnimationFrameId = 0;
let scrollSettleTimeoutId = 0;
let scrollLateTimeoutId = 0;
let glassRefreshNonce = 0;
let glassInitToken = 0;

const BASE_GLASS_CONFIG = {
  blurAmount: 0.5,
  refraction: 1.25,
  chromAberration: 0.25,
  edgeHighlight: 0.25,
  specular: 0,
  fresnel: 1,
  distortion: 0,
  cornerRadius: 18,
  zRadius: 18,
  saturation: 0,
  tintStrength: 0,
  brightness: -0.3,
  shadowSpread: 10,
  shadowOffsetY: 1
};

const MOBILE_MENU_GLASS_CONFIG = {
  ...BASE_GLASS_CONFIG,
  cornerRadius: 8,
  zRadius: 8,
  opacity: 1,
  shadowOpacity: 0.24
};

if (mobileMenuGlass) {
  mobileMenuGlass.dataset.config = JSON.stringify(MOBILE_MENU_GLASS_CONFIG);
}

document.querySelector("#current-year").textContent = new Date().getFullYear();

function updateHeaderGlass(forceRefresh = false) {
  const isScrolled = window.scrollY > 8;

  if (!forceRefresh && isScrolled === isHeaderScrolled && siteHeader.dataset.config) {
    return;
  }

  isHeaderScrolled = isScrolled;
  siteHeader.classList.toggle("is-scrolled", isScrolled);
  siteHeader.dataset.config = JSON.stringify({
    ...BASE_GLASS_CONFIG,
    opacity: isScrolled ? 1 : 0.02,
    edgeHighlight: isScrolled ? BASE_GLASS_CONFIG.edgeHighlight : 0,
    shadowOpacity: isScrolled ? 0.3 : 0,
    refreshNonce: forceRefresh ? glassRefreshNonce += 1 : glassRefreshNonce
  });
}

function handleScroll() {
  updateHeaderGlass(true);

  if (!scrollAnimationFrameId) {
    scrollAnimationFrameId = window.requestAnimationFrame(() => {
      liquidGlassInstance?.markChanged();
      scrollAnimationFrameId = 0;
    });
  }
  scheduleScrollGlassRefreshes();
}

function scheduleScrollGlassRefreshes() {
  window.clearTimeout(scrollSettleTimeoutId);
  window.clearTimeout(scrollLateTimeoutId);
  scrollSettleTimeoutId = window.setTimeout(refreshLiquidGlass, 90);
  scrollLateTimeoutId = window.setTimeout(refreshLiquidGlass, 240);
}

function refreshLiquidGlass() {
  window.requestAnimationFrame(() => {
    liquidGlassInstance?.markChanged(mainContent || undefined);
  });
}

function syncMobileMenuGlass() {
  if (!mobileMenuGlass || !navigationMenu.classList.contains("is-open")) {
    mobileMenuGlass?.classList.remove("is-open");
    return;
  }

  const rect = navigationMenu.getBoundingClientRect();
  mobileMenuGlass.classList.add("is-open");
  mobileMenuGlass.style.left = `${rect.left}px`;
  mobileMenuGlass.style.top = `${rect.top}px`;
  mobileMenuGlass.style.width = `${rect.width}px`;
  mobileMenuGlass.style.height = `${rect.height}px`;
}

function refreshMobileMenuGlass() {
  syncMobileMenuGlass();
  window.requestAnimationFrame(() => {
    liquidGlassInstance?.markChanged();
    window.setTimeout(() => liquidGlassInstance?.markChanged(), 80);
  });
}

function refreshGlassAfterAvatarReveal() {
  if (!profileAvatar?.src) {
    refreshLiquidGlass();
    return;
  }

  profileAvatar.decode()
    .catch(() => {})
    .finally(refreshLiquidGlass);
}

function syncThemeUI(theme) {
  const isDark = theme === "dark";
  const toggleLabel = isDark ? "Switch to light theme" : "Switch to dark theme";
  themeColorMeta.setAttribute("content", isDark ? "#202023" : "#ffffff");
  themeToggle.setAttribute("aria-label", toggleLabel);
  themeToggle.setAttribute("title", toggleLabel);
}

function setTheme(theme) {
  documentRoot.dataset.theme = theme;
  localStorage.setItem("theme", theme);
  syncThemeUI(theme);
  reinitializeLiquidGlass();
}

async function initLiquidGlass(initToken = ++glassInitToken) {
  if (!window.WebGLRenderingContext) {
    documentRoot.classList.add("liquid-glass-fallback");
    return;
  }

  try {
    await document.fonts.ready;
    const { LiquidGlass } = await import("https://cdn.jsdelivr.net/npm/@ybouane/liquidglass@1.0.3/dist/index.js");

    updateHeaderGlass();
    const nextGlassInstance = await LiquidGlass.init({
      root: glassRoot,
      glassElements: [siteHeader, mobileMenuGlass].filter(Boolean)
    });
    if (initToken !== glassInitToken) {
      nextGlassInstance.destroy();
      return;
    }
    liquidGlassInstance = nextGlassInstance;
    documentRoot.classList.add("liquid-glass-ready");
    syncMobileMenuGlass();
    refreshLiquidGlass();
  } catch {
    if (initToken === glassInitToken) {
      documentRoot.classList.add("liquid-glass-fallback");
    }
  }
}

function reinitializeLiquidGlass() {
  const initToken = ++glassInitToken;
  liquidGlassInstance?.destroy();
  liquidGlassInstance = undefined;
  documentRoot.classList.remove("liquid-glass-ready", "liquid-glass-fallback");
  updateHeaderGlass(true);
  void initLiquidGlass(initToken);
}

syncThemeUI(documentRoot.dataset.theme);
updateHeaderGlass();
void initLiquidGlass();

themeToggle.addEventListener("click", () => {
  setTheme(documentRoot.dataset.theme === "dark" ? "light" : "dark");
});

menuToggle.addEventListener("click", () => {
  const isOpen = menuToggle.getAttribute("aria-expanded") === "true";
  menuToggle.setAttribute("aria-expanded", String(!isOpen));
  menuToggle.setAttribute("aria-label", isOpen ? "Open navigation menu" : "Close navigation menu");
  navigationMenu.classList.toggle("is-open", !isOpen);
  refreshMobileMenuGlass();
});

navigationMenu.querySelectorAll("a").forEach((link) => {
  link.addEventListener("click", () => {
    navigationMenu.classList.remove("is-open");
    menuToggle.setAttribute("aria-expanded", "false");
    menuToggle.setAttribute("aria-label", "Open navigation menu");
    refreshMobileMenuGlass();
  });
});

identityContent?.addEventListener("pointerenter", () => {
  refreshGlassAfterAvatarReveal();
  window.setTimeout(refreshLiquidGlass, 260);
});

identityContent?.addEventListener("pointerleave", () => {
  refreshLiquidGlass();
  window.setTimeout(refreshLiquidGlass, 180);
});

window.addEventListener("resize", () => {
  if (window.innerWidth > 720) {
    navigationMenu.classList.remove("is-open");
    menuToggle.setAttribute("aria-expanded", "false");
  }
  refreshMobileMenuGlass();
  liquidGlassInstance?.markChanged();
});

window.addEventListener("scroll", handleScroll, { passive: true });

window.addEventListener("pagehide", () => {
  glassInitToken += 1;
  if (scrollAnimationFrameId) {
    window.cancelAnimationFrame(scrollAnimationFrameId);
  }
  window.clearTimeout(scrollSettleTimeoutId);
  window.clearTimeout(scrollLateTimeoutId);
  liquidGlassInstance?.destroy();
}, { once: true });
