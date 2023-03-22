function loading() {
    $("#loading").show();
    $("#content").hide();
}

window.onload = () => {
    const showSelectedDataTypeSection = (e) => {
        const radioButton = e.target;
        const selectedOption = radioButton.value;
        if (radioButton.checked) {
            switch (selectedOption) {
                case ("count-matrix"):
                    document.getElementById("data_type_warning").style.display = "none"
                    document.getElementById("count-matrix-input").style.display = "table-row";
                    document.getElementById("count_matrix").setAttribute("required", "")
                    document.getElementById("visium-10x-input").style.display = "none";
                    document.getElementById("visium_10x_matrix").removeAttribute("required")
                    document.getElementById("visium_10x_features").removeAttribute("required")
                    document.getElementById("visium_10x_barcodes").removeAttribute("required")
                    break;
                case ("visium-10x"):
                    document.getElementById("data_type_warning").style.display = "none"
                    document.getElementById("count-matrix-input").style.display = "none";
                    document.getElementById("count_matrix").removeAttribute("required")
                    document.getElementById("visium-10x-input").style.display = "table-row";
                    document.getElementById("visium_10x_matrix").setAttribute("required", "")
                    document.getElementById("visium_10x_features").setAttribute("required", "")
                    document.getElementById("visium_10x_barcodes").setAttribute("required", "")

                    break;
                default:
                    break;
            }
        }
    }

    const radioButtons = document.getElementById("data_types_selecting");

    radioButtons.addEventListener('change', showSelectedDataTypeSection);
}