/**
 * Design tokens for iLoveChem.
 *
 * Templates were written with stock Tailwind palettes (blue / indigo / purple /
 * sky / cyan / teal all mixed together). Instead of editing ~1,700 class names,
 * the palette is remapped here so every legacy class resolves to one accent
 * colour plus warm neutrals. Status colours (emerald / amber / red / rose) are
 * left alone because they carry meaning.
 */
const colors = require('tailwindcss/colors');

// The single accent. 600 is the interactive colour (white text passes AA).
const accent = {
  50: '#eff8f6', 100: '#d6eeea', 200: '#afddd6', 300: '#7ec6bc', 400: '#48a99d',
  500: '#238f82', 600: '#0f766e', 700: '#0d5f59', 800: '#0e4c48', 900: '#0f3f3b', 950: '#05211f',
};

module.exports = {
  content: ['../calculator/templates/**/*.html'],
  theme: {
    fontFamily: {
      sans: ['"IBM Plex Sans"', 'system-ui', '-apple-system', '"Segoe UI"', 'Roboto', 'sans-serif'],
      mono: ['"IBM Plex Mono"', 'ui-monospace', 'SFMono-Regular', 'Menlo', 'monospace'],
    },
    // Flat, consistent corners: one small radius for controls, one for panels.
    borderRadius: {
      none: '0', sm: '3px', DEFAULT: '4px', md: '6px', lg: '6px',
      xl: '8px', '2xl': '8px', '3xl': '8px', full: '9999px',
    },
    // Borders do the separating. Only overlays get a real shadow.
    boxShadow: {
      none: 'none',
      '2xs': 'none', xs: 'none', sm: 'none', DEFAULT: 'none', inner: 'none',
      md: '0 2px 8px rgb(20 20 18 / 0.08)',
      lg: '0 4px 16px rgb(20 20 18 / 0.10)',
      xl: '0 8px 24px rgb(20 20 18 / 0.12)',
      '2xl': '0 12px 32px rgb(20 20 18 / 0.16)',
    },
    // "backdrop-blur-xs" is a Tailwind v4 name that does not exist in v3.
    backdropBlur: { none: '0', xs: '0', sm: '0', DEFAULT: '0', md: '0' },
    extend: {
      colors: {
        slate: colors.stone,
        gray: colors.stone,
        zinc: colors.stone,
        neutral: colors.stone,
        blue: accent,
        indigo: accent,
        violet: accent,
        purple: accent,
        sky: accent,
        cyan: accent,
        teal: accent,
        fuchsia: accent,
        pink: colors.rose,
        paper: '#f7f6f3',
        ink: '#1c1b19',
      },
      fontSize: {
        // Nothing smaller than 12px: the old [10px] micro-labels were unreadable.
        '2xs': ['0.75rem', { lineHeight: '1rem' }],
      },
    },
  },
  plugins: [],
};
