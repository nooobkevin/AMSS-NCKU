這是AMSS-NCKU程序：
```
AMSS - NCKU is a numerical relativity program developed in China, which is used to numerically solve Einstein's equations and calculate the change of the gravitational field over time.

AMSS - NCKU uses the finite difference method and the adaptive mesh refinement technique to achieve the numerical solution of Einstein's equations.

Currently, AMSS - NCKU can successfully handle binary black hole systems and multiple black hole systems, calculate the time evolution of these systems, and solve the gravitational waves released during these processes.
```
其中amss-ncku-python/AMSS_NCKU_source是原代碼（2007年寫的），其他的應該是interface。已知這是一坨屎山。根據以下建議，搞懂這個repo到底要計算什麼（核心方程），以什麼方法（方程動力學），然後提供重寫建議（比如，已知裡面有求解Einstein Field Equation的非常重複的代碼，可以用python的einsum解決）。
```
第一階段：搞懂「它想幹什麼」（輸入與接口）
不要先看 C/Fortran 代碼，先看 Python 外殼，弄清楚這個程序需要什麼參數。

amss-ncku-python/AMSS_NCKU_Input.py (以及 inputfile_example/ 下的文件)

目的：這是控制台。看這裡你可以知道程序接受哪些物理參數（黑洞質量、自旋、距離）和數值參數（網格大小、時間步長）。
關注點：變量名。這些變量名通常會直接對應到 C++ 內部的變量，是你搜索代碼的關鍵詞。
amss-ncku-python/AMSS_NCKU_Program.py

目的：了解整個工作流。
關注點：它如何調用編譯好的二進制文件？傳遞了什麼命令行參數？
第二階段：尋找入口與主循環（骨架）
進入 AMSS_NCKU_source 目錄。你需要找到 main() 函數在哪裡。

setup.C (極大可能是入口)

推測：在沒有明顯 main.C 的情況下，setup.C 通常包含程序的初始化邏輯，甚至包含 main 函數本身。
關注點：尋找 int main(...)。看它是如何讀取參數，以及如何初始化 cgh (Computational Grid Hierarchy) 的。
driver.h / bssn_class.C (演化核心)

背景：BSSN (Baumgarte-Shapiro-Shibata-Nakamura) 是數值相對論中最常用的方程形式。
關注點：
在 bssn_class.C 中尋找類似 evolve(), step(), 或 run() 的函數。這是時間演化的心臟。
看它是如何調用 Runge-Kutta 積分器 (rungekutta4_rout.f90) 的。
第三階段：物理模塊（血肉）
現在你可以看具體的物理實現了。

TwoPunctures.C / TwoPunctures.h (初始條件)

背景：這是雙黑洞模擬最經典的初始數據方法（Puncture method）。
關注點：如何設置初始的度規（Metric）和曲率（Extrinsic Curvature）。這是模擬的起點。
bssn_rhs.f90 (右端項方程)

警告：這文件裡全是數學公式的 Fortran 實現，非常枯燥且容易看不懂。
關注點：不要一行行讀。只需要知道這裡是計算愛因斯坦方程的 $\partial_t \phi = \dots$ 部分即可。它是被 C++ 的 bssn_class 調用的。
第四階段：基礎設施（地基 - 最容易暈的地方）
除非你要修 Bug 或改並行邏輯，否則初期盡量少看這部分，因為這裡通常是「屎山」最臭的地方（手寫的內存管理、MPI 通信、網格拼接）。

cgh.C / cgh.h

這應該是 Computational Grid Hierarchy。管理整個網格結構。
patch.C / MPatch.C / ShellPatch.C

AMR（自適應網格細化）的核心。處理不同分辨率網格之間的插值。AMSS-NCKU 似乎使用了多個 Patch 拼接的技術（尤其是 ShellPatch 可能用於處理波區）。
Parallel.C / Parallel_bam.C

MPI 並行通信的實現。
總結推薦閱讀路徑
熱身：inputfile_example/BSSN_Input.py (看懂輸入參數)
入口：grep "main" *.C (找到主函數，很可能是 setup.C)
核心對象：bssn_class.h (看這個類定義了哪些成員變量，理解數據結構)
時間循環：bssn_class.C (看 step 函數如何推進時間)
物理方程：bssn_rhs.f90 (瀏覽一下即可，確認是算 BSSN 方程)
初始數據：TwoPunctures.C
避坑指南
忽略 .cu 文件：那是 GPU 加速代碼，邏輯和 CPU 版一樣但更難讀，先看 CPU 版。
忽略 Z4c 系列：Z4c 是另一種方程形式（比 BSSN 新），先看懂 BSSN，Z4c 的結構是一樣的。
忽略 Ansorg.C / ABE.C：這些可能是特定的求解器或舊代碼，先抓主幹。
```

過程中，你可以利用DRAFTSHEET.md來放一些你的草稿。每次有進展，都把你對這個repo的認知更新到NOTE.md。