// Keep the theme's private jQuery instance, as expected by bootstrap-sphinx.js.
window.$jqTheme = jQuery.noConflict(true);

// Recent Sphinx versions no longer provide a global jQuery. Bootstrap 3's
// version check still reads that global, even though its plugins use $jqTheme.
if (typeof window.jQuery === "undefined") {
  window.jQuery = window.$jqTheme;
}
