example = 0

N = 3
h = [5.7, 5.0][example]  # mm
w = [52, 50][example] # mm
crack_width = 34.29787 - 16.516718 # mm, narrow version
# crack_width = 34.534954 - 16.358664 # mm, wide version
load = 100 # kg
# fs = 1.823 # kg/mm из статьи Л. В. Степановой и В. С. Долгих
fs = 1.945659 # kg/mm из программы тарировки
# fs = 1.77

sigma_22_inf = load / (w*h) # kg/mm**2

# K, M = 5, 3
K, M = 7, 5
# K, M = 15, 0
# K, M = 15, 9

path = f"data/{'One_crack' if example==0 else 'Two cracks'}/"

image_file = [
    f"{path}Р01-{load}кг.jpg",
    f"{path}Р 11-{load}кг.jpg"
    ][example]

points_file = [
    f"{path}points_Р01-{load}кг_N_{N}_of_5.txt",
    f"{path}Points_Р 11-{load}кг_N_{N}_of_4.txt"
    ][example]

center = [
    [34.29787, 51.604862],
    [38.214905, 31.455807]
    ][example]

# center = [34.3769, 51.604862]
