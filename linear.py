def LinRegress(*args):
    xlist = args[0]
    ylist = args[1]
    InitialGuess = args[2]
    if not isinstance(xlist, (np.ndarray, np.generic)):
        xlist = np.array(xlist)
    if not isinstance(ylist, (np.ndarray, np.generic)):
        ylist = np.array(ylist)
    LineFunc = lambda k,x,m: k * x + m
    Error = 1
    NewError = 0.5
    i = 0
    while Error > NewError:
        try:
            Slopes = [((ylist[i+1] - ylist[i])/(xlist[i+1] - xlist[i])) for i in range(0, len(ylist)-1)]
            Errors = [Slopes[i]* i - ylist[i] for i in range(0, len(ylist)-1)]
            print(Slopes)


        except Exception as E:
            raise E
        i += 1
        if i > 20:
            print("Did not converge within given limit")
            return None
