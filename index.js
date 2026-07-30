const root = document.documentElement;
const glassRoot = document.querySelector("#glass-root");
const header = document.querySelector(".site-header");
const contentRoot = document.querySelector("main");
const themeToggle = document.querySelector(".theme-toggle");
const menuToggle = document.querySelector(".menu-toggle");
const menu = document.querySelector(".nav-links");
const themeMeta = document.querySelector('meta[name="theme-color"]');
const identityCopy = document.querySelector(".identity-copy");
const profileAvatar = document.querySelector("[data-profile-avatar]");

let liquidGlassInstance;
let headerIsScrolled = false;
let scrollFrame = 0;
let scrollSettleTimer = 0;
let scrollLateTimer = 0;
let liquidRefreshNonce = 0;
let liquidGlassInitToken = 0;
let liquidGlassInitPromise;

const GLASS_CONFIG = {
  blurAmount: 0.25,
  refraction: 1.25,
  chromAberration: 0.05,
  edgeHighlight: 0.05,
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

document.querySelector("#current-year").textContent = new Date().getFullYear();

function updateHeaderGlass(forceRefresh = false) {
  const isScrolled = window.scrollY > 8;

  if (!forceRefresh && isScrolled === headerIsScrolled && header.dataset.config) {
    return;
  }

  headerIsScrolled = isScrolled;
  header.classList.toggle("is-scrolled", isScrolled);
  header.dataset.config = JSON.stringify({
    ...GLASS_CONFIG,
    opacity: isScrolled ? 1 : 0.02,
    edgeHighlight: isScrolled ? GLASS_CONFIG.edgeHighlight : 0,
    shadowOpacity: isScrolled ? 0.3 : 0,
    refreshNonce: forceRefresh ? liquidRefreshNonce += 1 : liquidRefreshNonce
  });
}

function handleScroll() {
  updateHeaderGlass(true);

  if (scrollFrame) {
    window.clearTimeout(scrollSettleTimer);
    window.clearTimeout(scrollLateTimer);
    scrollSettleTimer = window.setTimeout(refreshLiquidGlass, 90);
    scrollLateTimer = window.setTimeout(refreshLiquidGlass, 240);
    return;
  }

  scrollFrame = window.requestAnimationFrame(() => {
    liquidGlassInstance?.markChanged();
    scrollFrame = 0;
  });
  window.clearTimeout(scrollSettleTimer);
  window.clearTimeout(scrollLateTimer);
  scrollSettleTimer = window.setTimeout(refreshLiquidGlass, 90);
  scrollLateTimer = window.setTimeout(refreshLiquidGlass, 240);
}

function refreshLiquidGlass() {
  window.requestAnimationFrame(() => {
    liquidGlassInstance?.markChanged(contentRoot || undefined);
  });
}

function refreshLiquidGlassAfterAvatarLoad() {
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
  themeMeta.setAttribute("content", isDark ? "#202023" : "#ffffff");
  themeToggle.setAttribute("aria-label", toggleLabel);
  themeToggle.setAttribute("title", toggleLabel);
}

function setTheme(theme) {
  root.dataset.theme = theme;
  localStorage.setItem("theme", theme);
  syncThemeUI(theme);
  reinitializeLiquidGlass();
}

async function initLiquidGlass(token = ++liquidGlassInitToken) {
  if (!window.WebGLRenderingContext) {
    root.classList.add("liquid-glass-fallback");
    return;
  }

  try {
    await document.fonts.ready;
    const { LiquidGlass } = await import("https://cdn.jsdelivr.net/npm/@ybouane/liquidglass@1.0.3/dist/index.js");

    updateHeaderGlass();
    const nextInstance = await LiquidGlass.init({
      root: glassRoot,
      glassElements: [header]
    });
    if (token !== liquidGlassInitToken) {
      nextInstance.destroy();
      return;
    }
    liquidGlassInstance = nextInstance;
    root.classList.add("liquid-glass-ready");
    refreshLiquidGlass();
  } catch {
    if (token === liquidGlassInitToken) {
      root.classList.add("liquid-glass-fallback");
    }
  }
}

function reinitializeLiquidGlass() {
  const token = ++liquidGlassInitToken;
  liquidGlassInstance?.destroy();
  liquidGlassInstance = undefined;
  root.classList.remove("liquid-glass-ready", "liquid-glass-fallback");
  updateHeaderGlass(true);
  liquidGlassInitPromise = initLiquidGlass(token);
}

syncThemeUI(root.dataset.theme);
updateHeaderGlass();
liquidGlassInitPromise = initLiquidGlass();

themeToggle.addEventListener("click", () => {
  setTheme(root.dataset.theme === "dark" ? "light" : "dark");
});

menuToggle.addEventListener("click", () => {
  const isOpen = menuToggle.getAttribute("aria-expanded") === "true";
  menuToggle.setAttribute("aria-expanded", String(!isOpen));
  menuToggle.setAttribute("aria-label", isOpen ? "Open navigation menu" : "Close navigation menu");
  menu.classList.toggle("is-open", !isOpen);
});

menu.querySelectorAll("a").forEach((link) => {
  link.addEventListener("click", () => {
    menu.classList.remove("is-open");
    menuToggle.setAttribute("aria-expanded", "false");
    menuToggle.setAttribute("aria-label", "Open navigation menu");
  });
});

identityCopy?.addEventListener("pointerenter", () => {
  refreshLiquidGlassAfterAvatarLoad();
  window.setTimeout(refreshLiquidGlass, 260);
});

identityCopy?.addEventListener("pointerleave", () => {
  refreshLiquidGlass();
  window.setTimeout(refreshLiquidGlass, 180);
});

window.addEventListener("resize", () => {
  if (window.innerWidth > 720) {
    menu.classList.remove("is-open");
    menuToggle.setAttribute("aria-expanded", "false");
  }
  liquidGlassInstance?.markChanged();
});

window.addEventListener("scroll", handleScroll, { passive: true });

window.addEventListener("pagehide", () => {
  liquidGlassInitToken += 1;
  if (scrollFrame) {
    window.cancelAnimationFrame(scrollFrame);
  }
  window.clearTimeout(scrollSettleTimer);
  window.clearTimeout(scrollLateTimer);
  liquidGlassInstance?.destroy();
}, { once: true });
