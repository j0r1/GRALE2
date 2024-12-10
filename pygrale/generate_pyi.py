import sys

def _commentsremoved(l):
    newL = []
    for x in l:
        if x == "#":
            break
        newL.append(x)
    return "".join(newL)

def _preprocess(lines):
    newLines = []
    linesToMerge = None
    parenthCount = 0
    for l in lines:
        if l.strip().startswith("def ") or l.strip().startswith("cdef "):
            assert parenthCount == 0, "Expecting parenthCount to be zero"

            for x in l:
                if x == '(':
                    parenthCount += 1
                elif x == ')':
                    parenthCount -= 1
                    assert parenthCount >= 0, "Negative parenthesis count!"

            if parenthCount == 0:
                newLines.append(l)
            else:
                assert linesToMerge is None
                linesToMerge = [ _commentsremoved(l).rstrip() ]

        else:
            if parenthCount == 0: # Not in a 'def' expression
                newLines.append(l)
                assert linesToMerge is None
            else:
                assert linesToMerge

                for x in l:
                    if x == '(':
                        parenthCount += 1
                    elif x == ')':
                        parenthCount -= 1
                        assert parenthCount >= 0, "Negative parenthesis count!"
        
                if parenthCount == 0: # final parenthesis closed
                    linesToMerge.append(" " + _commentsremoved(l.strip()))
                    newLines.append(" ".join(linesToMerge) + "\n")
                    linesToMerge = None
                else: # still not closed
                    linesToMerge.append(" " + _commentsremoved(l.strip()))


    #print("".join(newLines))
    return newLines

def splitIntoArguments(argsStr):
    parentCounts = { "()": 0, "[]": 0 }
    arg = ""

    for idx in range(len(argsStr)):
        c = argsStr[idx]
        if c == "(":
            parentCounts["()"] += 1
            arg += c
        elif c == ")":
            parentCounts["()"] -= 1
            arg += c
        elif c == "[":
            parentCounts["[]"] += 1
            arg += c
        elif c == "]":
            parentCounts["[]"] -= 1
            arg += c
        elif c == ",":
            if parentCounts["()"] == 0 and parentCounts["[]"] == 0:
                return [ arg.strip() ] + splitIntoArguments(argsStr[idx+1:])
            
            arg += c
        else:
            arg += c

    return [ arg.strip() ]


def gen_pyi(inFn, outFn):
    if type(inFn) == list:
        lines = []
        for fn in inFn:
            lines += open(fn).readlines()
    else:
        lines = open(inFn).readlines()

    lines = _preprocess(lines)

    state = {
        "needStart": True,
        "secondPrevLine": None,
        "prevLine": None,
        "commentLines": [],
        "codeComments": [],
        "modComments": None,
        "lastIndent": 0
    }

    def findNotSpace(l):
        idx = 0
        while idx < len(l) and l[idx] == " ":
            idx += 1
        return idx
    
    def hasStart(l):
        return l.strip().startswith('"""') or l.strip().startswith('r"""')
    
    def hasStartAndEnd(l: str):
        idx1 = l.find('"""')
        assert idx1 >= 0, "Internal error: couldn't find start anymore"
        idx2 = l.find('"""', idx1+3)
        return True if idx2 >= 0 else False
    
    def changeCDef(l:str):
        if l.strip().startswith("cdef class "):
            idx = l.find("cdef class ")
            return l[:idx] + "class " + l[idx+11:]

        if l.strip().startswith("cdef "):
            idx = l.find("cdef ")
            return l[:idx] + "def " + l[idx+5:]
        
        return l

    def removeSquareBrackets(l):
        newLine = ""
        prevPos = 0
        while True:
            idx1 = l.find("[", prevPos)
            if idx1 < 0:
                return newLine + l[prevPos:]
            
            idx2 = l.find("]", idx1)
            assert idx2 > idx1, "Need closing ]"

            newLine += l[prevPos:idx1]
            prevPos = idx2+1

    def filterTypes(arg, typeNames):
        for t in typeNames:
            if arg.startswith(t):
                return arg[len(t):]
        return arg

    def changeArgs(l:str):
        idx1 = l.find("(")
        if idx1 < 0:
            return l
        
        startPart = l[:idx1+1]
        idx2 = l.rfind(")")
        assert idx2 > idx1, "Parentheses not closed"
        endPart = l[idx2:]

        middle = l[idx1+1:idx2]
        #print("TODO: consider arguments:")
        #print(splitIntoArguments(middle))

        #middle = removeSquareBrackets(middle)
        
        #parts = middle.split(",")
        parts = splitIntoArguments(middle)
        newParts = []
        for arg in parts:
            arg = arg.strip()
            arg = filterTypes(arg, ["cbool ", "double ", "np.ndarray ", "bool ", "int ",
                                    "bytes ", "LensPlane ", "np.ndarray[np.float32_t,ndim=1] ",
                                    "np.ndarray[double, ndim=2] ", "np.ndarray[double,ndim=2] "])
            newParts.append(arg)

        middle = ", ".join(newParts)

        return startPart + middle + endPart
        

    def processCommentLines():
        if state["prevLine"] == None:
            assert state["modComments"] is None, "Internal error"
            state["modComments"] = state["commentLines"]
        else:
            if state["secondPrevLine"].strip().startswith("@"):
                state["codeComments"].append(state["secondPrevLine"])
            state["codeComments"].append(changeArgs(changeCDef(state["prevLine"])))
            for c in state["commentLines"]:
                state["codeComments"].append(c)
            state["codeComments"].append(" "*(state["lastIndent"] + 4) + "...\n")
            state["codeComments"].append("\n")

        state["commentLines"] = []
     
    for l in lines:
        if state["needStart"]:
            if hasStart(l):
                state["commentLines"].append(l)
                if hasStartAndEnd(l):
                    processCommentLines()
                else:
                    state["needStart"] = False
            else:
                state["secondPrevLine"] = state["prevLine"]
                state["prevLine"] = l
                state["lastIndent"] = findNotSpace(l)
        else: # not needStart, in comment
            state["commentLines"].append(l)
            if '"""' in l:
                state["needStart"] = True
                processCommentLines()

    if state["modComments"] is None:
        state["modComments"] = []

    open(outFn,"wt").write(''.join(state["modComments"] + [ "\n" ] + state["codeComments"]))

def main():
    gen_pyi(sys.argv[1:-1], sys.argv[-1])

if __name__ == "__main__":
    main()
