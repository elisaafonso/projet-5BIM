#make movie
import subprocess
import shutil 
import os 

if __name__ == "__main__": 
    path_results = "/home/vidium06/src/ANA_SENS/Results_PhysiCell/sensitivity_analysis_descent_time2"
    video_path = os.path.normpath(os.path.join(path_results, "output_video"))
    os.makedirs(video_path, exist_ok=True)
    
    for i in range (64) : 
        path_output_i = os.path.join(path_results, f"output_{i}")
        process1 = subprocess.run(
            [
                "make",
                "jpeg",
                f"OUTPUT={path_output_i}"
            ],
            capture_output=True,
            text=True,
            cwd="/home/vidium06/src/ANA_SENS/PhysiCell"
        )

        # Affichier les outputs et erreurs
        print("Output:")
        print(process1.stdout) 
        print("Error:")
        print(process1.stderr)
        
        process2 = subprocess.run(
            [
                "make",
                "movie",
                f"OUTPUT={path_output_i}"
            ],
            capture_output=True,
            text=True,
            cwd="/home/vidium06/src/ANA_SENS/PhysiCell"
        )

        # Affichier les outputs et erreurs
        print("Output:")
        print(process2.stdout) 
        print("Error:")
        print(process2.stderr)
        shutil.move(os.path.join(path_results, f"output_{i}", "out.mp4"), os.path.join(video_path, f"video_{i}"))

