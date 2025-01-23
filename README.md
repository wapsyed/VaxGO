# VaxGO: An interactive application for systems vaccinology analysis of transcriptomic data

**Authors:**
Wasim Aluísio Prates-Syed1,2, Aline Aparecida Lima¹, Nelson Cortes¹, Evelyn Carvalho¹, Jaqueline Dinis Queiroz Silva¹, Bárbara Hamaguchi¹, José Krieger4, Thomas Hagan³, Gustavo Cabral-Miranda¹.
**Institutions:**
1. Biomedical Sciences Institute, University of São Paulo, São Paulo, Brazil.
2. Pro-Vaccine Union, University of São Paulo, Ribeirão Preto, São Paulo, Brazil.
3. Cincinnati Children's Hospital, Cincinnati, Ohio, USA.
4. Heart Institute (INCOR), University of São Paulo.

**About VaxGO**
RNAseq is crucial for investigating transcriptional patterns in immunology, especially within immunization research. However, RNAseq data analysis often requires programming skills, which can restrict access for researchers without such expertise.We introduce VaxGO, a tool designed to facilitate the analysis of differentially expressed genes in the context of immune processes and cells during vaccination. This tool integrates data from Gene Ontology, CellMarker 2.0, the MSigDB Vax collection, and other key studies, including transcriptional atlases of vaccines against COVID-19 and other diseases. VaxGO is an interactive, web-based tool, offering a user-friendly platform for exploring immune responses and vaccine efficacy without programming expertise.


(This tool is still in development. It may be slow on the browser, so we recommend downloading this repository and running the "VaxGO_tool.rmd" on RStudio, as described in more details below).
**Access the tool here:** https://wapsyed.shinyapps.io/VaxGO_Tool_test

**How to Download and Set Up the Tool:**
1) **Install R and RStudio**
If you don’t already have them, download and install R and the Desktop version of RStudio.

2)**Clone the Repository**
- Open RStudio.
- On the top-left menu, click File → New Project.
- Select Version Control, then Git.
- Paste the repository URL: https://github.com/wapsyed/VaxGO.git and click Create Project. The project will open automatically after the download is complete.

3) **Pre-Install Dependencies**
- Locate the file preinstall.r in the project folder.
- Open the file in RStudio and click Run to execute the script. This will install all the necessary packages.

4) **Run the VaxGO Tool**
- Open the file VaxGO_Tool.rmd.
- Click the Run Document button in RStudio to launch the tool.
