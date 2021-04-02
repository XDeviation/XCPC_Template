
import os
path = os.listdir(os.getcwd())

used = []
with open('complete.md', 'r') as f:
    for line in f.readlines():
        used.append(line.strip())

for p in path:
    if os.path.isdir(p):
        if '.' in p:
            continue
        if p in used:
            continue
        print(p)
        