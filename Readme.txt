---------------------------------------------------------------------------

OCT Denoising package for volumetric OCT images

---------------------------------------------------------------------------
 Contents
---------------------------------------------------------------------------

The package contains these files

*) demo_cirrus3D_deSpeckle.m        : denoising demo script
*) demo_cirrus3D_deSpeckle_GPU.m    : denoising using GPU acceleration 
                                        demo script
*) demo_batch_process.m             : batch demo
*) deSpeckle_batch_process.m        : a batch process function that denoise
                                        all the .img files in the folder
*) install_path.m                   : add all the required files in Matlab
                                        searching path

---------------------------------------------------------------------------
 Installation & Usage
---------------------------------------------------------------------------

Unzip OCT_DeSpeckle_Package.zip (contains codes) in a folder that is in the 
MATLAB path. Execute the script "install_path.m" to install all the required 
packages. 

Execute the script "demo_cirrus3D_deSpeckle.m" to run the denoise demo, or 
execute the script "demo_cirrus3D_deSpeckle.m" to run the GPU demo.
Execute the script "demo_batch_process.m' ro run a batch demo that process 
all the .img files in specific folder automatically.

You can freely modify the parameters involved in the filtering at the 
beginning of each demo.

---------------------------------------------------------------------------
 Requirements
---------------------------------------------------------------------------

*) MS Windows 64 bit, Linux 64 bit or Mac OS X 64 bit
*) Matlab R2011b or later with installed:
   -- Image Processing Toolbox 
   -- Signal Processing Toolbox 
   -- Wavelet Toolbox 
   -- Statistics and Machine Learning Toolbox (Only for automatically 
        parameter setting)
   -- Parallel Computing Toolbox (Only for GPU acceleration)

-------------------------------------------------------------------