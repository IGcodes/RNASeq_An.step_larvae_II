document.addEventListener("DOMContentLoaded", function() {

    // --- MANUALLY CONFIGURE YOUR FILES HERE ---
    const fileData = {
        "GO_enrichment": {
            displayName: "GO Enrichment",
            files: [
                { displayName: "Adama vs Erer (Up Regulated)", filePath: "data/GO_enrichment/sigUp_Adama_vs_Erer_GO.csv", type: "GO" },
                { displayName: "Adama vs Erer (Down Regulated)", filePath: "data/GO_enrichment/sigDown_Adama_vs_Erer_GO.csv", type: "GO" },
                { displayName: "Jijiga vs Erer (Up Regulated)", filePath: "data/GO_enrichment/sigDown_Erer_vs_Jijiga_GO.csv", type: "GO" },
                { displayName: "Jijiga vs Erer (Down Regulated)", filePath: "data/GO_enrichment/sigUp_Erer_vs_Jijiga_GO.csv", type: "GO" },
                { displayName: "Adama vs Jijiga (Up Regulated)", filePath: "data/GO_enrichment/sigUp_Adama_vs_Jijiga_GO.csv", type: "GO" },
                { displayName: "Adama vs Jijiga (Down Regulated)", filePath: "data/GO_enrichment/sigDown_Adama_vs_Jijiga_GO.csv", type: "GO" },
                { displayName: "Wild vs Lab (Up Regulated)", filePath: "data/GO_enrichment/sigUp_Wild_vs_Lab_GO.csv", type: "GO" },
                { displayName: "Wild vs Lab (Down Regulated)", filePath: "data/GO_enrichment/sigDown_Wild_vs_Lab_GO.csv", type: "GO" }
            ]
        },
        "KEGGS": {
            displayName: "KEGG Pathways",
            files: [
                // { displayName: "Adama vs Erer (Down Regulated)", filePath: "data/KEGGS/Down_Adama_vs_Erer_KEGG.csv", type: "KEGG" }
                { displayName: "Adama vs Erer (Up Regulated)", filePath: "data/KEGGS/Up_Adama_vs_Erer_KEGG.csv", type: "KEGG" },
                { displayName: "Adama vs Erer (Down Regulated)", filePath: "data/KEGGS/Down_Adama_vs_Erer_KEGG.csv", type: "KEGG" },
                { displayName: "Jijiga vs Erer (Up Regulated)", filePath: "data/KEGGS/Down_Erer_vs_Jijiga_KEGG.csv", type: "KEGG" },
                { displayName: "Jijiga vs Erer (Down Regulated)", filePath: "data/KEGGS/Up_Erer_vs_Jijiga_KEGG.csv", type: "KEGG" },
                { displayName: "Adama vs Jijiga (Up Regulated)", filePath: "data/KEGGS/Up_Adama_vs_Jijiga_KEGG.csv", type: "KEGG" },
                { displayName: "Adama vs Jijiga (Down Regulated)", filePath: "data/KEGGS/Down_Adama_vs_Jijiga_KEGG.csv", type: "KEGG" },
                { displayName: "Wild vs Lab (Up Regulated)", filePath: "data/KEGGS/Up_Wild_vs_Lab_KEGG.csv", type: "KEGG" },
                { displayName: "Wild vs Lab (Down Regulated)", filePath: "data/KEGGS/Down_Wild_vs_Lab_KEGG.csv", type: "KEGG" }
            ]
        },
        "GSEA": {
            displayName: "Gene Set Enrichment Analysis",
            files: [
                { displayName: "Adama vs Erer Results", filePath: "data/GSEA/Adama_vs_Erer_GSEA_results.csv", type: "GSEA" },
                { displayName: "Adama vs Jijiga Results", filePath: "data/GSEA/Adama_vs_Jijiga_GSEA_results.csv", type: "GSEA" },
                { displayName: "Erer vs Jijiga Results", filePath: "data/GSEA/Erer_vs_Jijiga_GSEA_results.csv", type: "GSEA" },
                { displayName: "Wild vs Lab Results", filePath: "data/GSEA/Wild_vs_Lab_GSEA_results.csv", type: "GSEA" }
            ]
        }
    };

    // --- Get references to the HTML elements ---
    const directorySelect = document.getElementById("directory-select");
    const fileSelect = document.getElementById("file-select");
    const tableContainer = document.getElementById("table-container");
    const chartCanvas = document.getElementById("bar-chart");
    let chartInstance = null; // To hold the chart object

    // --- CSV Parsing Function ---
    function parseCSV(text) {
        const lines = text.trim().split('\n');
        const header = lines[0].split(',').map(h => h.trim().replace(/"/g, ''));
        const rows = lines.slice(1).map(line => {
            const values = line.split(',').map(v => v.trim().replace(/"/g, ''));
            let obj = {};
            header.forEach((h, i) => {
                obj[h] = values[i];
            });
            return obj;
        });
        return rows;
    }

    // --- Table Creation Function ---
    function createTable(data) {
        if (!data || data.length === 0) {
            tableContainer.innerHTML = "<p>No data to display.</p>";
            return;
        }
        const headers = Object.keys(data[0]);
        let tableHTML = '<table><thead><tr>';
        headers.forEach(h => tableHTML += `<th>${h}</th>`);
        tableHTML += '</tr></thead><tbody>';
        data.forEach(row => {
            tableHTML += '<tr>';
            headers.forEach(h => tableHTML += `<td>${row[h]}</td>`);
            tableHTML += '</tr>';
        });
        tableHTML += '</tbody></table>';
        tableContainer.innerHTML = tableHTML;
    }

    // --- Chart Creation Function ---
    function createChart(data, fileType) {
        if (chartInstance) {
            chartInstance.destroy(); // Destroy the old chart before creating a new one
        }
        const ctx = chartCanvas.getContext('2d');
        let chartConfig = {
            type: 'bar',
            data: {
                labels: [],
                datasets: [{
                    label: '',
                    data: [],
                    backgroundColor: 'rgba(54, 162, 235, 0.6)',
                    borderColor: 'rgba(54, 162, 235, 1)',
                    borderWidth: 1
                }]
            },
            options: {
                indexAxis: 'y', // To make bars horizontal
                responsive: true,
                maintainAspectRatio: false,
                scales: { y: { beginAtZero: true } },
                plugins: { legend: { display: false }, title: { display: true, text: '' } }
            }
        };

        // Customize chart based on file type
        switch (fileType) {
            case 'GO':
            case 'KEGG':
                // Plot Top 15 by 'Count'
                const sortedGOData = data.sort((a, b) => b.Count - a.Count).slice(0, 15).reverse();
                chartConfig.data.labels = sortedGOData.map(row => row.Description);
                chartConfig.data.datasets[0].data = sortedGOData.map(row => row.Count);
                chartConfig.data.datasets[0].label = 'Gene Count';
                chartConfig.options.plugins.title.text = 'Top 15 Enriched Terms by Gene Count';
                break;

            case 'GSEA':
                // Plot Top 15 by absolute 'NES'
                const sortedGSEAData = data.sort((a, b) => Math.abs(b.NES) - Math.abs(a.NES)).slice(0, 15).reverse();
                chartConfig.data.labels = sortedGSEAData.map(row => row.Description);
                chartConfig.data.datasets[0].data = sortedGSEAData.map(row => row.NES);
                chartConfig.data.datasets[0].label = 'Normalized Enrichment Score (NES)';
                chartConfig.options.plugins.title.text = 'Top 15 GSEA Results by Normalized Enrichment Score';
                // Color bars based on positive/negative NES
                chartConfig.data.datasets[0].backgroundColor = sortedGSEAData.map(row => row.NES > 0 ? 'rgba(255, 99, 132, 0.6)' : 'rgba(54, 162, 235, 0.6)');
                chartConfig.data.datasets[0].borderColor = sortedGSEAData.map(row => row.NES > 0 ? 'rgba(255, 99, 132, 1)' : 'rgba(54, 162, 235, 1)');
                break;
        }

        chartInstance = new Chart(ctx, chartConfig);
    }

    // --- Main Function to Fetch and Display Data ---
    async function displayData() {
        const selectedOption = fileSelect.options[fileSelect.selectedIndex];
        const filePath = selectedOption.value;
        const fileType = selectedOption.dataset.type;

        if (!filePath) {
            tableContainer.innerHTML = '<p>Select a file to see its data here.</p>';
            if (chartInstance) chartInstance.destroy();
            return;
        }

        try {
            const response = await fetch(filePath);
            if (!response.ok) throw new Error(`HTTP error! Status: ${response.status}`);
            const text = await response.text();
            
            const parsedData = parseCSV(text);
            createTable(parsedData);
            createChart(parsedData, fileType);

        } catch (error) {
            console.error("Error loading file:", error);
            tableContainer.innerHTML = `<p>Error loading file: ${filePath}. See console for details.</p>`;
            if (chartInstance) chartInstance.destroy();
        }
    }

    // --- Functions to manage dropdowns ---
    function populateDirectorySelect() {
        directorySelect.innerHTML = `<option value="">-- Select Analysis Type --</option>`;
        for (const key in fileData) {
            const option = document.createElement("option");
            option.value = key;
            option.textContent = fileData[key].displayName;
            directorySelect.appendChild(option);
        }
    }

    function updateFileSelect() {
        const selectedDirKey = directorySelect.value;
        fileSelect.innerHTML = `<option value="">-- Select a file --</option>`;
        if (selectedDirKey && fileData[selectedDirKey]) {
            fileData[selectedDirKey].files.forEach(file => {
                const option = document.createElement("option");
                option.value = file.filePath;
                option.textContent = file.displayName;
                option.dataset.type = file.type; // Store the file type in the option
                fileSelect.appendChild(option);
            });
        }
        // Clear displays
        tableContainer.innerHTML = '<p>Select a file to see its data here.</p>';
        if (chartInstance) chartInstance.destroy();
    }

    // --- Event Listeners ---
    directorySelect.addEventListener("change", updateFileSelect);
    fileSelect.addEventListener("change", displayData);

    // --- Initial setup ---
    populateDirectorySelect();
});