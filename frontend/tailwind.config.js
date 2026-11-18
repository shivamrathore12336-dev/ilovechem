/**
 * Professional Scientific Design Tokens for iLoveChem (Swiss Lab Workbench Edition).
 */
const colors = require('tailwindcss/colors');

const accent = {
  50:  '#effaf7',
  100: '#d7f3ec',
  200: '#b0e6da',
  300: '#7cd2c2',
  400: '#48b7a5',
  500: '#249c8b',
  600: '#0f766e',
  700: '#0d5f59',
  800: '#0e4c48',
  900: '#0f3f3b',
  950: '#052321',
};

module.exports = {
  content: ['../calculator/templates/**/*.html'],
  theme: {
    fontFamily: {
      display: ['"Space Grotesk"', '"IBM Plex Sans"', 'system-ui', 'sans-serif'],
      sans: ['"IBM Plex Sans"', 'system-ui', '-apple-system', '"Segoe UI"', 'Roboto', 'sans-serif'],
      mono: ['"IBM Plex Mono"', 'ui-monospace', 'SFMono-Regular', 'Menlo', 'monospace'],
    },
    borderRadius: {
      none: '0',
      sm: '4px',
      DEFAULT: '6px',
      md: '8px',
      lg: '10px',
      xl: '12px',
      '2xl': '14px',
      '3xl': '16px',
      full: '9999px',
    },
    boxShadow: {
      none: 'none',
      '2xs': '0 1px 2px 0 rgb(15 23 42 / 0.03)',
      xs: '0 1px 2px 0 rgb(15 23 42 / 0.04)',
      sm: '0 1px 3px 0 rgb(15 23 42 / 0.06), 0 1px 2px -1px rgb(15 23 42 / 0.04)',
      DEFAULT: '0 2px 6px -1px rgb(15 23 42 / 0.07), 0 1px 3px -1px rgb(15 23 42 / 0.05)',
      md: '0 4px 12px -2px rgb(15 23 42 / 0.08), 0 2px 6px -2px rgb(15 23 42 / 0.04)',
      lg: '0 10px 24px -4px rgb(15 23 42 / 0.10), 0 4px 10px -3px rgb(15 23 42 / 0.05)',
      xl: '0 16px 32px -6px rgb(15 23 42 / 0.12), 0 6px 14px -4px rgb(15 23 42 / 0.06)',
      '2xl': '0 24px 48px -12px rgb(15 23 42 / 0.18)',
      inner: 'inset 0 2px 4px 0 rgb(15 23 42 / 0.04)',
    },
    extend: {
      colors: {
        slate: colors.slate,
        gray: colors.slate,
        zinc: colors.slate,
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
        paper: '#f6f8fa',
        ink: '#0f172a',
      },
      fontSize: {
        '2xs': ['0.75rem', { lineHeight: '1.1rem' }],
      },
    },
  },
  plugins: [],
};
