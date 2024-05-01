import pickle

scores = {}
with open("src/experiment/data/batch-stress-scores.pkl", "rb") as f:
    scores = pickle.load(f)

for score in scores.items():
    print(score)