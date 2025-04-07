import numpy as np
import matplotlib.pyplot as plt
import pickle as pl
import h5py

def writeplot(fig_handle,file):
    print(f'Number of axes in the figure: {len(fig_handle.axes)}')
    ax=fig_handle.axes[0]
    print(f'Number of lines in the axis: {len(ax.lines)}')

    list = []
    for i in range(len(ax.lines)):
        line = ax.lines[i] 
        x=np.array(line.get_xdata())
        y=np.array(line.get_ydata())
        st=line.get_linestyle()
        co=line.get_color()
        data = np.column_stack((x,y))
        list.append((data,st,co))

#    with open(file+'.pkl', 'wb') as f:
#        pl.dump(list, f)
#        print(f'Written to file:',file+'.pkl')

    with h5py.File(file+'.h5', 'w') as f:
        f['header'] = "Matplotlib figure data V01"
        f['axis000/count'] = len(list)
        f['axis000/lim'] = (ax.get_xlim(),ax.get_ylim())
        f['axis000/scale'] = (ax.get_xscale(),ax.get_yscale())
        f['axis000/label'] = (ax.get_xlabel(),ax.get_ylabel())
        
        for i in range(len(list)):
            nn = 'axis000/line'+str(i).zfill(3)
            f.create_dataset(nn+'/data',data=list[i][0])
            f.create_dataset(nn+'/st',data=list[i][1])
            f.create_dataset(nn+'/co',data=list[i][2])
        print(f'Written to file:',file+'.h5')

def readplot_pkl(file):
    with open(file+'.pkl', 'rb') as f:
        loaded_list = pl.load(f)
        print(f'Read file:',file,' with ',len(loaded_list),' lines.')
    return loaded_list

def readplot(file):
    with h5py.File(file+'.h5', 'r') as f:
        header = f['header'][()]
        header=header.decode('utf-8')

        lim   = f['axis000/lim'][()]
        scale = f['axis000/scale'][()]
        label = f['axis000/label'][()]
        list1 = (header,lim,scale,label)

        ll = f['axis000/count'][()]
        list = []
        for i in range(ll):
            nn = 'axis000/line'+str(i).zfill(3)
            da = f[nn+'/data'][:]
            st = f[nn+'/st'][()]
            st=st.decode('utf-8')
            co = f[nn+'/co'][()]
            co=co.decode('utf-8')
            list.append((da,st,co))
        print('Loaded '+header + ": "+str(ll)+" lines")
    return list,list1

def showplot(file,show=True):
    d,ll=readplot(file)
    fig_handle=plt.figure()
    for i in range(0,len(d)):
        data = d[i][0]
        st = d[i][1]
        co = d[i][2]
        plt.plot(data[:,0],data[:,1],linestyle=st,color=co)
    plt.xscale(ll[2][0].decode('utf-8'))
    plt.yscale(ll[2][1].decode('utf-8'))
    plt.xlabel(ll[3][0].decode('utf-8'))
    plt.ylabel(ll[3][1].decode('utf-8'))
    plt.xlim(ll[1][0]);plt.ylim(ll[1][1])
    if show:
        plt.show()
    return fig_handle