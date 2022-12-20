window.onload = function() {
    const fill_vals = () => {
        document.getElementById('num_std_thr').value = '0.5';
        document.getElementById('mu').value = '0.95';
        document.getElementById('path_len').value = '3';
        document.getElementById('num_iter').value = '10000';
    }
    document.getElementById('fill_default_vals').addEventListener('click', fill_vals)
}