// Change these three values when you want to customize the site's identity.
const SITE_CONFIG = {
  name: "KWEIWEN TSENG",
  title: "KWEIWEN TSENG",
};

const root = document.documentElement;
const glassRoot = document.querySelector("#glass-root");
const header = document.querySelector(".site-header");
const themeToggle = document.querySelector(".theme-toggle");
const menuToggle = document.querySelector(".menu-toggle");
const menu = document.querySelector(".nav-links");
const themeMeta = document.querySelector('meta[name="theme-color"]');
const identityCopy = document.querySelector(".identity-copy");

let liquidGlassInstance;
let headerIsScrolled = false;
let scrollFrame = 0;

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

document.querySelectorAll("[data-site-name]").forEach((element) => {
  element.textContent = SITE_CONFIG.name;
});
const heroTitle = document.querySelector("[data-hero-title]");
const heroSubtitle = document.querySelector("[data-hero-subtitle]");

if (heroTitle) {
  heroTitle.textContent = SITE_CONFIG.title;
}
if (heroSubtitle) {
  heroSubtitle.textContent = SITE_CONFIG.subtitle;
}
document.querySelector("#current-year").textContent = new Date().getFullYear();

function updateHeaderGlass() {
  const isScrolled = window.scrollY > 8;

  if (isScrolled === headerIsScrolled && header.dataset.config) {
    return;
  }

  headerIsScrolled = isScrolled;
  header.classList.toggle("is-scrolled", isScrolled);
  header.dataset.config = JSON.stringify({
    ...GLASS_CONFIG,
    opacity: isScrolled ? 1 : 0.02,
    edgeHighlight: isScrolled ? GLASS_CONFIG.edgeHighlight : 0,
    shadowOpacity: isScrolled ? 0.3 : 0
  });
}

function handleScroll() {
  updateHeaderGlass();

  if (scrollFrame) {
    return;
  }

  scrollFrame = window.requestAnimationFrame(() => {
    liquidGlassInstance?.markChanged();
    scrollFrame = 0;
  });
}

function refreshLiquidGlass() {
  window.requestAnimationFrame(() => {
    liquidGlassInstance?.markChanged();
  });
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
  liquidGlassInstance?.markChanged();
}

async function initLiquidGlass() {
  if (!window.WebGLRenderingContext) {
    root.classList.add("liquid-glass-fallback");
    return;
  }

  try {
    await document.fonts.ready;
    const { LiquidGlass } = await import("https://cdn.jsdelivr.net/npm/@ybouane/liquidglass@1.0.3/dist/index.js");

    updateHeaderGlass();
    liquidGlassInstance = await LiquidGlass.init({
      root: glassRoot,
      glassElements: [header]
    });
    root.classList.add("liquid-glass-ready");
  } catch {
    root.classList.add("liquid-glass-fallback");
  }
}

syncThemeUI(root.dataset.theme);
updateHeaderGlass();
initLiquidGlass();

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
  refreshLiquidGlass();
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
  if (scrollFrame) {
    window.cancelAnimationFrame(scrollFrame);
  }
  liquidGlassInstance?.destroy();
}, { once: true });
