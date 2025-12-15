# DMM: Diffraction Matrix Method

Matlab code implementation of the Diffraction Matrix Method for calculating the optical properties of periodically corrugated OLEDs.

* **Authors:** Chanhyung Park, Jeongmin Shin, Min Seok Jang
* **Last Update:** 2025.12.15

---

## 📌 Notice (Dependencies)

This code requires the following external tools and toolboxes to function correctly and efficiently.

### 1. RETICOLO (RCWA Solver)
This code utilizes **RETICOLO**, a rigorous coupled-wave analysis (RCWA) open-source code developed by J. P. Hugonin and P. Lalanne.
* **Required Action:**
  1. Please download **RETICOLO V8** or a higher version.
  2. Unzip the RETICOLO folder to the same directory where `DMM_main.m` or `DMM_2D_main.m` is located.

### 2. MATLAB Parallel Computing Toolbox (⚡ Critical for Speed)
This code is optimized using **`parfor`** (parallel for-loops) to handle heavy matrix computations.

> [!IMPORTANT]
> **⚠️ PERFORMANCE WARNING:**
> It is **strongly recommended** to install and use the **MATLAB Parallel Computing Toolbox**.
>
> * **With Toolbox:** The code runs in parallel, fully utilizing all your CPU cores.
> * **Without Toolbox:** The code defaults to serial processing, which can be **30 to 50 times slower**.
>
> **Please ensure this toolbox is installed to avoid extremely long calculation times.**
---

## ▶️ How to Run This Code

Run the corresponding main file in MATLAB:

* For **1D corrugated OLEDs**:
    ```matlab
    DMM_main.m
    ```
* For **2D corrugated OLEDs**:
    ```matlab
    DMM_2D_main.m
    ```

---

## Citation

> Chanhyung Park, Jeongmin Shin, Sanmun Kim, Songju Lee, Juho Park, Jaehyeok Park, Sehong Park, Seunghyup Yoo, and Min Seok Jang, "Fast and rigorous optical simulation of periodically corrugated light-emitting diodes based on a diffraction matrix method," Opt. Express 31, 20410-20423 (2023)
>
> **Link:** https://opg.optica.org/oe/fulltext.cfm?uri=oe-31-12-20410
