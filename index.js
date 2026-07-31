const root = document.documentElement;
const glassRoot = document.querySelector("#glass-root");
const header = document.querySelector(".site-header");
const mobileMenuGlass = document.querySelector(".mobile-menu-glass");
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

const MENU_GLASS_CONFIG = {
  ...GLASS_CONFIG,
  cornerRadius: 8,
  zRadius: 8,
  opacity: 1,
  shadowOpacity: 0.24
};

if (mobileMenuGlass) {
  mobileMenuGlass.dataset.config = JSON.stringify(MENU_GLASS_CONFIG);
}

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

function syncMobileMenuGlass() {
  if (!mobileMenuGlass || !menu.classList.contains("is-open")) {
    mobileMenuGlass?.classList.remove("is-open");
    return;
  }

  const rect = menu.getBoundingClientRect();
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
      glassElements: [header, mobileMenuGlass].filter(Boolean)
    });
    if (token !== liquidGlassInitToken) {
      nextInstance.destroy();
      return;
    }
    liquidGlassInstance = nextInstance;
    root.classList.add("liquid-glass-ready");
    syncMobileMenuGlass();
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
  refreshMobileMenuGlass();
});

menu.querySelectorAll("a").forEach((link) => {
  link.addEventListener("click", () => {
    menu.classList.remove("is-open");
    menuToggle.setAttribute("aria-expanded", "false");
    menuToggle.setAttribute("aria-label", "Open navigation menu");
    refreshMobileMenuGlass();
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
  refreshMobileMenuGlass();
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
