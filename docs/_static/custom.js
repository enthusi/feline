insertCurrentYearToPageFooter = () => {
    document.addEventListener("DOMContentLoaded", () => {
    let contentInfo = document.querySelector('[role="contentinfo"] p');
    if (contentInfo) {
        contentInfo.innerHTML = "&#169; Copyright " + new Date().getFullYear() + ", Feline Team.";
    }
});
}

insertCurrentYearToPageFooter();
