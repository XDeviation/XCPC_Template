
import os
path = os.listdir(os.getcwd())
for p in path:
    if os.path.isdir(p):
        pathdir = "./" + p
        newpath = os.listdir(pathdir)
        for fil in newpath:
            if 'cpp' in fil:
                print(' ')
                print('```cpp')
                print('//'+ pathdir + '/' + fil)
                os.system(f'cat -b {pathdir}/{fil} > tmp.md')
                fo = open(f"tmp.md", "r+")
                tmpstr = fo.read()
                print(tmpstr)
                fo.close()
                os.system(f'rm -rf tmp.md')
                print(' ')
                print('```')
        