import os

def list_python_and_config_files(folder="."):
    print("=== PROJECT PYTHON & SCRIPT FILES ===")
    ignored_dirs = {"venv", ".venv", ".git", "__pycache__", "maf_files", 
                    "chl_sccc_2023", "coad_silu_2022", "lihc_tcga_pan_can_atlas_2018"}
    
    for root, dirs, files in os.walk(folder):
        dirs[:] = [d for d in dirs if d not in ignored_dirs]
        for file in files:
            if file.endswith(('.py', '.sh', '.bat', '.json', '.yaml', '.yml')) and file != "diagnose.py":
                rel_path = os.path.relpath(os.path.join(root, file), folder)
                print(f"📄 {rel_path}")

if __name__ == "__main__":
    list_python_and_config_files()