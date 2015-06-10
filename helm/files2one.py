import os

path = r'C:\Users\yinz\Desktop\Penlight-master'
last_pos = path.rfind('\\')+1

fp = open('all.txt', 'w')
for root, dirs, files in os.walk(path):
    for f in files:
        #if f.startswith('.'): continue
        text = ''.join(open(os.path.join(root,f),'r').readlines())
        header = '####$$$$ %s/%s'%(root[last_pos:], f)
        fp.write(header)
        fp.write('\n')
        fp.write(text)
        
fp.close()
