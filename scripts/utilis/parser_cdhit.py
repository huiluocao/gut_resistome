sample_host={}
with open('./1174_sample.txt','r') as f:
    for line in f.readlines():
        sample_host.update({line.split('\t')[0]:line.split('\t')[2]})
sample_host

with open('./arg_1174_ant6_800_cdhit.fasta.clstr','r') as f:
    chunks = re.split(r'\n(?=>\b *)', f.read())
    for chunk in chunks:
        hosts=[]
        lines=chunk.split('\n')
        for line in lines:
            if '*' in line:
                #seq=(line.split('$')[1])[:-5]
                seq=(line.split('>')[1])[:-5]
                sample=(line.split('$'))[0].split('>')[1]
                host=sample_host[sample]
                #print(seq,sample)
                hosts.append(host)
            elif "Cluster" in line:
                pass
            else:
                try:
                    sample=(line.split('$')[0]).split('>')[1]
                    host=sample_host[sample]
                    hosts.append(host)                    
                except:
                    samples=''
        #print(seq,'\t'.join(samples))
        with open('arg_1174_ant6_800_cdhit_host_count.txt','a') as g:
            g.write(seq+'\t'+str(hosts.count('Human'))+'\t'+str(hosts.count('Swine'))+'\n')
