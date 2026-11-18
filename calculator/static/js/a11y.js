/*
 * Shared accessibility, CSRF cookie retrieval, 3Dmol WebGL diagnostics,
 * and RFC 4180 Blob download helpers for iLoveChem.
 */
(function () {
    function getCookie(name) {
        if (!document.cookie) return null;
        var cookies = document.cookie.split(';');
        for (var i = 0; i < cookies.length; i++) {
            var cookie = cookies[i].trim();
            if (cookie.substring(0, name.length + 1) === (name + '=')) {
                return decodeURIComponent(cookie.substring(name.length + 1));
            }
        }
        return null;
    }

    function getCsrfToken() {
        var fromCookie = getCookie('csrftoken');
        if (fromCookie) return fromCookie;
        var input = document.querySelector('input[name="csrfmiddlewaretoken"]');
        return input ? input.value : '';
    }

    function hasWebGLSupport() {
        try {
            var canvas = document.createElement('canvas');
            return !!(
                window.WebGLRenderingContext &&
                (canvas.getContext('webgl') || canvas.getContext('experimental-webgl'))
            );
        } catch (e) {
            return false;
        }
    }

    function renderWebGLFallback(container, message) {
        if (!container) return;
        var emptyOverlay = container.querySelector('#viewer-empty-state') ||
                           document.getElementById('viewer-empty-state');
        if (emptyOverlay) emptyOverlay.classList.add('hidden');

        var fallback = container.querySelector('.webgl-fallback-notice');
        if (!fallback) {
            fallback = document.createElement('div');
            fallback.className = 'webgl-fallback-notice absolute inset-0 flex flex-col items-center justify-center p-6 text-center bg-slate-900/95 text-slate-200 z-20';
            container.appendChild(fallback);
        }
        fallback.innerHTML =
            '<div class="max-w-sm space-y-2.5">' +
            '<div class="inline-flex items-center px-2.5 py-1 rounded text-xs font-mono bg-amber-500/20 text-amber-300 border border-amber-500/30">3D WebGL Unavailable</div>' +
            '<p class="text-sm text-slate-300 leading-relaxed">' + (message || '3D coordinates were generated cleanly, but hardware WebGL rendering is unavailable in this browser. Use the 2D Structure view or download the 3D coordinate file directly.') + '</p>' +
            '</div>';
    }

    function downloadBlob(content, filename, mimeType) {
        var blob = new Blob([content], { type: mimeType || 'text/plain;charset=utf-8;' });
        var url = URL.createObjectURL(blob);
        var link = document.createElement('a');
        link.href = url;
        link.setAttribute('download', filename);
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
        setTimeout(function () { URL.revokeObjectURL(url); }, 250);
    }

    function enhance(root) {
        (root || document).querySelectorAll('div[onclick], span[onclick], li[onclick]').forEach(function (el) {
            if (el.hasAttribute('tabindex')) return;
            el.setAttribute('role', 'button');
            el.setAttribute('tabindex', '0');
            el.addEventListener('keydown', function (e) {
                if (e.key === 'Enter' || e.key === ' ') {
                    e.preventDefault();
                    el.click();
                }
            });
        });
    }

    function init() {
        enhance(document);
        if (window.lucide) lucide.createIcons();
    }

    window.getCookie = getCookie;
    window.getCsrfToken = getCsrfToken;
    window.enhanceClickables = enhance;
    window.ilovechemUtils = {
        getCookie: getCookie,
        getCsrfToken: getCsrfToken,
        hasWebGLSupport: hasWebGLSupport,
        renderWebGLFallback: renderWebGLFallback,
        downloadBlob: downloadBlob
    };

    if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
    else init();
})();
