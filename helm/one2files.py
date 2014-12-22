import os

one = r'd:\home\tools\all.txt'
target = r'd:\home\tools'

def get_filename(fn):
    subs = fn.split('/')
    curr = target
    for sub in subs[0:-1]:
        curr = os.path.join(curr, sub)
        if not os.path.exists(curr):
            os.makedirs(curr)
    return os.path.join(curr, subs[-1])
            

text = ''.join(open(one,'r').readlines())
texts = text.split('####$$$$')
for text in texts:
    text = text.strip()
    if len(text) == 0: continue
    lines = text.split('\n')
    filename = lines[0]
    filename = filename.replace('\\', '/')
    filename = get_filename(filename)
    with open(filename,'w') as f:
        f.write('\n'.join(lines[1:]))

