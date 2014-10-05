from config import PICKLE_PROTOCOL
import pickle

def f_save(obj, file_name):
    w_stream = open(file_name, 'wb')
    pickle.dump(obj, w_stream, PICKLE_PROTOCOL)
    w_stream.close()
    return "\nContent saved to file " + file_name

def f_recover(file_name):
    r_stream = open(file_name,'rb')
    thawed = pickle.load(r_stream)
    r_stream.close()
    return thawed