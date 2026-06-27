(() => {
  const scriptUrl = (() => {
    try {
      const current = document.currentScript && document.currentScript.src
        ? document.currentScript.src
        : "pwa-enhancements.js";
      return new URL(current, window.location.href);
    } catch {
      return new URL(window.location.href);
    }
  })();

  const serviceWorkerUrl = new URL("./sw.js", scriptUrl);
  const appIconUrl = new URL("./icons/icon-192.png", scriptUrl).href;
  const installHintDismissKey = "wms_rot_install_hint_dismissed";

  const isAppleDevice = () => {
    const platform = navigator.platform || "";
    const userAgent = navigator.userAgent || "";
    const touchPoints = Number(navigator.maxTouchPoints || 0);
    return /iPhone|iPad|iPod|Mac/i.test(platform)
      || /iPhone|iPad|iPod|Macintosh/i.test(userAgent)
      || (platform === "MacIntel" && touchPoints > 1);
  };

  const isStandaloneMode = () => {
    const displayModeStandalone = window.matchMedia
      && window.matchMedia("(display-mode: standalone)").matches;
    return Boolean(displayModeStandalone || navigator.standalone === true);
  };

  const isLikelySafari = () => {
    const userAgent = navigator.userAgent || "";
    const vendor = navigator.vendor || "";
    if (!/Apple/i.test(vendor)) return false;
    return !/CriOS|FxiOS|EdgiOS|OPiOS|DuckDuckGo|YaBrowser|SamsungBrowser|Chrome|Chromium/i.test(userAgent);
  };

  const supportsWebAppInstall = () => {
    if (!window.isSecureContext) return false;
    if (!("serviceWorker" in navigator)) return false;

    // Chromium-style install prompt support.
    if ("BeforeInstallPromptEvent" in window) return true;

    // Apple Safari install flow (Share -> Add to Home Screen / Add to Dock).
    if (isAppleDevice() && isLikelySafari()) return true;

    return false;
  };

  let swRegistrationPromise = null;
  const registerServiceWorker = () => {
    if (swRegistrationPromise) return swRegistrationPromise;
    swRegistrationPromise = (async () => {
      if (!("serviceWorker" in navigator) || !window.isSecureContext) return null;
      try {
        const scope = new URL("./", serviceWorkerUrl).pathname;
        return await navigator.serviceWorker.register(serviceWorkerUrl.href, { scope });
      } catch (error) {
        console.debug("WMS-Rot service worker registration skipped:", error);
        return null;
      }
    })();
    return swRegistrationPromise;
  };

  let storagePersistencePromise = null;
  const requestPersistentStorage = () => {
    if (storagePersistencePromise) return storagePersistencePromise;
    storagePersistencePromise = (async () => {
      if (!navigator.storage || typeof navigator.storage.persist !== "function") return false;
      try {
        const isAlreadyPersisted = typeof navigator.storage.persisted === "function"
          ? await navigator.storage.persisted()
          : false;
        if (isAlreadyPersisted) return true;
        return await navigator.storage.persist();
      } catch (error) {
        console.debug("WMS-Rot storage persistence request skipped:", error);
        return false;
      }
    })();
    return storagePersistencePromise;
  };

  let notificationPermissionPromise = null;
  const prepareNotificationsFromGesture = async () => {
    void requestPersistentStorage();
    if (!("Notification" in window)) return false;
    if (Notification.permission === "granted") return true;
    if (Notification.permission === "denied") return false;
    if (!notificationPermissionPromise) {
      notificationPermissionPromise = Notification.requestPermission()
        .then((permission) => permission === "granted")
        .catch(() => false)
        .finally(() => {
          notificationPermissionPromise = null;
        });
    }
    return notificationPermissionPromise;
  };

  const shouldNotifyNow = (force) => {
    if (force) return true;
    const hasFocusFn = typeof document.hasFocus === "function";
    return document.hidden || (hasFocusFn && !document.hasFocus());
  };

  const notifyBackground = async ({
    title = document.title || "WMS-Rot",
    body = "",
    tag = "wms-rot",
    force = false,
  } = {}) => {
    if (!shouldNotifyNow(force)) return false;
    if (!("Notification" in window)) return false;
    if (Notification.permission !== "granted") return false;

    const options = {
      body: String(body || ""),
      tag: String(tag || "wms-rot"),
      renotify: true,
      icon: appIconUrl,
      badge: appIconUrl,
      data: { url: window.location.href },
      timestamp: Date.now(),
    };

    try {
      const registration = await registerServiceWorker();
      if (registration && typeof registration.showNotification === "function") {
        await registration.showNotification(String(title || "WMS-Rot"), options);
      } else {
        new Notification(String(title || "WMS-Rot"), options);
      }
      return true;
    } catch (error) {
      console.debug("WMS-Rot notification skipped:", error);
      return false;
    }
  };

  const showAppleInstallHint = () => {
    if (!supportsWebAppInstall()) return;
    if (!isAppleDevice() || isStandaloneMode()) return;
    try {
      if (window.localStorage.getItem(installHintDismissKey) === "1") return;
    } catch {
      // no-op
    }

    const hint = document.createElement("div");
    hint.setAttribute("role", "status");
    hint.style.position = "fixed";
    hint.style.left = "12px";
    hint.style.right = "12px";
    hint.style.bottom = "calc(12px + env(safe-area-inset-bottom, 0px))";
    hint.style.zIndex = "1200";
    hint.style.display = "flex";
    hint.style.alignItems = "center";
    hint.style.gap = "10px";
    hint.style.padding = "10px 12px";
    hint.style.borderRadius = "10px";
    hint.style.border = "1px solid rgba(15, 23, 42, 0.2)";
    hint.style.background = "rgba(15, 23, 42, 0.95)";
    hint.style.color = "#e2e8f0";
    hint.style.fontSize = "0.92rem";
    hint.style.boxShadow = "0 8px 30px rgba(2, 6, 23, 0.35)";

    const message = document.createElement("span");
    message.textContent = "Apple/iOS: open Share and choose Add to Home Screen to install the full-screen web app.";
    message.style.flex = "1 1 auto";

    const closeButton = document.createElement("button");
    closeButton.type = "button";
    closeButton.textContent = "Chiudi";
    closeButton.style.border = "1px solid rgba(148, 163, 184, 0.4)";
    closeButton.style.background = "transparent";
    closeButton.style.color = "inherit";
    closeButton.style.borderRadius = "8px";
    closeButton.style.padding = "4px 8px";
    closeButton.style.cursor = "pointer";
    closeButton.addEventListener("click", () => {
      hint.remove();
      try {
        window.localStorage.setItem(installHintDismissKey, "1");
      } catch {
        // no-op
      }
    });

    hint.appendChild(message);
    hint.appendChild(closeButton);
    document.body.appendChild(hint);
  };

  window.WMSPwaEnhancements = Object.freeze({
    isAppleDevice,
    isStandaloneMode,
    supportsWebAppInstall,
    registerServiceWorker,
    requestPersistentStorage,
    prepareNotificationsFromGesture,
    notifyBackground,
  });

  void registerServiceWorker();
  void requestPersistentStorage();
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", showAppleInstallHint, { once: true });
  } else {
    showAppleInstallHint();
  }
})();
