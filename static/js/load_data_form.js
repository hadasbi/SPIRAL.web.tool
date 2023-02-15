function loading() {
    $("#loading").show();
    $("#content").hide();
}

window.onload = () => {
    const showSelectedDataTypeSection = (e) => {
        console.log('hi! in here with button', e.target)
        const radioButton = e.target;
        const selectedOption = radioButton.value;
        console.log('selecet', selectedOption)
        if (radioButton.checked) {
            switch (selectedOption) {
                case ("count-matrix"):
                    document.getElementById("data_type_warning").style.display = "none"
                    document.getElementById("count-matrix-input").style.display = "table-row";
                    document.getElementById("visium-10x-input").style.display = "none";
                    break;
                case ("visium-10x"):
                    document.getElementById("data_type_warning").style.display = "none"
                    document.getElementById("count-matrix-input").style.display = "none";
                    document.getElementById("visium-10x-input").style.display = "table-row";
                    break;
                default:
                    break;
            }
        }
    }

    const radioButtons = document.getElementById("data_types_selecting");

    radioButtons.addEventListener('change', showSelectedDataTypeSection);
}