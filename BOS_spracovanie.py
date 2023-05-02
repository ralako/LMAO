import BOS

# Hmotnost jednej castice (kg)
mC_start = 1e-18
mC_stop = 1e-12
mC_num = 50


# Hmotnost oblaku (g)
m_start = 1
m_stop = 1000
m_num = 50


# Vzdialenost od stred oblaku (cm)
Dy_start = 0
Dy_stop = 100
Dy_num = 50


# Gulicka s priemerom 1cm
mD = 1.36    #g
C = 0.5     #
S = 0.0000785  #m^2
v = 13000   #m/s


#-----Spracovanie-----#
odpad = BOS.Odpad(mD,C,S,v)
#odpad.prelet_stredom(mC_start,mC_stop,mC_num, m_start,m_stop,m_num)
odpad.prelet_bokom(Dy_start,Dy_stop,Dy_num, m_start, m_stop,m_num,mC=1e-16)
