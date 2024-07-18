import sys

def gen_pyi(inFn, outFn):
    if type(inFn) == list:
        lines = []
        for fn in inFn:
            lines += open(fn).readlines()
    else:
        lines = open(inFn).readlines()

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
        middle = removeSquareBrackets(middle)
        
        parts = middle.split(",")
        newParts = []
        for arg in parts:
            arg = arg.strip()
            arg = filterTypes(arg, ["cbool ", "double ", "np.ndarray ", "bool ", "int ",
                                    "bytes ", "LensPlane "])
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
