# VaxGO: An interactive application for systems vaccinology analysis of transcriptomic data

**Authors:**

Wasim Aluísio Prates-Syed1,2, Aline Aparecida Lima¹, Nelson Cortes¹, Evelyn Carvalho¹, Jaqueline Dinis Queiroz Silva¹, Bárbara Hamaguchi¹, José Krieger4, Thomas Hagan³, Gustavo Cabral-Miranda¹.
**Institutions:**

1. Biomedical Sciences Institute, University of São Paulo, São Paulo, Brazil.
2. Pro-Vaccine Union, University of São Paulo, Ribeirão Preto, São Paulo, Brazil.
3. Cincinnati Children's Hospital, Cincinnati, Ohio, USA.
4. Heart Institute (INCOR), University of São Paulo.

**About VaxGO**

RNA sequencing (RNAseq) is crucial for investigating transcriptional patterns in immunology and vaccine research. However, the analysis of RNAseq data often requires programming skills, which can limit accessibility for researchers lacking such expertise. We present VaxGO, an intuitive web-based tool designed to facilitate the analysis of differentially expressed genes in the context of immune processes and cells during vaccination. This tool integrates data from Gene Ontology, CellMarker 2.0, the MSigDB Vax collection, and other key studies, including transcriptional atlases of vaccines against COVID-19 and other diseases. VaxGO is an interactive, web-based tool, offering a user-friendly platform for exploring immune responses and vaccine efficacy without programming expertise. 

**How to use this tool**

**Browser version:** [Click here](https://wapsyed.shinyapps.io/VaxGO_Shiny/).

This tool is still in development. It may be slow on the browser, so we recommend downloading this repository and running it locally with the "VaxGO_tool.rmd" file on RStudio, as described in more details below.

**How to download and use it on RStudio:**

1) **Install R and RStudio**
If you don’t already have them, download and install R and the Desktop version of RStudio.

2) **Clone the repository**
- Open RStudio.
- On the top-left menu, click File → New Project.
- Select Version Control, then Git.
- Paste the repository URL: https://github.com/wapsyed/VaxGO.git and click Create Project. The project will open automatically after the download is complete.

3) **Pre-install dependencies**
- Locate the file "preinstall.r" in the project folder.
- Open the file in RStudio and click Run to execute the script. This will install all the necessary packages.

4) **Run the VaxGO tool**
- Open the file "VaxGO_Tool.rmd".
- Click the "Run Document" button in RStudio to launch the tool.

5) **Run an example**
- Click on the button "Upload" and upload the file "Example_DEGs_BNT_multiple.csv" from the "Example" folder. This will open all DEGs from conditions of the vaccination with BNT162b2/Comirnaty described in "COVID-19 Vaccination Atlas". 
