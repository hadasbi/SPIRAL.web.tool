window.onload = function() {
    const fill_vals = () => {
        //document.getElementById('num_std_thr').value = '0.5';
        var values="0.5,0.75,1";
        $.each(values.split(","), function(i,e){
            $("#num_std_thr option[value='" + e + "']").prop("selected", true);
        });
        document.getElementById('mu').value = '0.95';
        document.getElementById('path_len').value = '3';
        document.getElementById('num_iter').value = '10000';
    }
    document.getElementById('fill_default_vals').addEventListener('click', fill_vals)
}