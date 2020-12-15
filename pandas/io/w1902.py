import sys

# map01.py
for line in sys.stdin:
    data = line.split(",")
    if len(data) == 3:
        print(data[0] + "\t" + "{},{}".format(data[1], data[2]))
    else:
        print(data[1] + "\t" + "{},{},{}".format(data[2], data[3], 1))


# reduce01.py
weight = {}
votes = {}
for line in sys.stdin:
    data = line.split("\t")
    key = data[0]
    value = data[1].split(",")
    if len(value) == 2:
        if key not in weight:
            weight[key] = {}
        weight[key][value[1]] = value[0]

    else:
        if key not in votes:
            votes[key] = {value[1]: {value[0]: 1}}
        else:
            votes[key][value[1]][value[0]] += 1

for state in votes.keys():
    for course in votes[state].keys():
        for food in votes[state][course].keys():
            print(
                "{}\t{},{},{}".format(
                    state,
                    course,
                    food,
                    votes[state][course][food] * weight[state][course],
                )
            )

# map02.py
for line in sys.stdin:
    data = line.split("\t")
    key = data[0]
    value = data[1].split(",")
    print("{}\t{},{}".format(value[0], value[1], value[2]))


# reduce02.py
course_food = {}
for line in sys.stdin:
    data = line.split("\t")
    course = data[0]
    value = data[1].split(",")
    food = value[0]
    weighted_sum = value[1]
    if course not in course_food:
        course_food[course] = {"food": food, "weighted_sum": weighted_sum}
    else:
        if weighted_sum > course_food[course]["weighted_sum"]:
            course_food[course]["food"] = food
            course_food[course]["weighted_sum"] = weighted_sum


for course in course_food.keys():
    print("{},{}".format(course, course_food[course]["food"]))
