import os

def ensurePathExists(folder):
    if not os.path.exists(folder):
        os.makedirs(folder)

module_path = os.path.dirname(os.path.abspath(__file__))

par_storage = os.path.relpath(os.path.join(module_path, "cache"), ".") # .par file storage path
tim_storage = os.path.relpath(os.path.join(module_path, "cache"), ".") # .tim file storage path
out_storage = os.path.relpath(os.path.join(module_path, "cache"), ".") # .out file storage path

ensurePathExists(par_storage)
ensurePathExists(tim_storage)
ensurePathExists(out_storage)

if __name__ == "__main__":
    raise Exception("This is a settings file, no point in running it.")
