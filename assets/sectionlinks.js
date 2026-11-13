// Documenter renders a sidebar section ("Applications") as a header that only opens and closes the
// section; make its name a link to the section's first page (its overview). The chevron still toggles.
document.addEventListener("DOMContentLoaded", function () {
  document.querySelectorAll(".docs-menu .tocitem > .docs-label, .docs-menu span.tocitem").forEach(function (label) {
    var header = label.classList.contains("tocitem") ? label : label.parentElement;
    var list = header.nextElementSibling;
    var first = list && list.querySelector("a.tocitem");
    if (!first) return;
    var a = document.createElement("a");
    a.href = first.getAttribute("href");
    a.textContent = label.textContent;
    a.style.color = "inherit";
    a.addEventListener("click", function (e) { e.stopPropagation(); });   // navigate, don't toggle
    label.textContent = "";
    label.appendChild(a);
  });
});
