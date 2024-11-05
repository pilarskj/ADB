package adbp;

// class to store all information about a branch
public class Branch {
    public int startNode;
    public int endNode;
    public double startTime;
    public double endTime;
    public String branchMode;
    public int branchIndex;
    public int leftIndex;
    public int rightIndex;
    public int nodeType; // type at branch end

    public Branch(int startNode, int endNode, double startTime, double endTime, String branchMode) {
        // assertions
        assert branchMode.equals("internal") || branchMode.equals("external");
        assert startTime < endTime;

        this.startNode = startNode;
        this.endNode = endNode;
        this.startTime = startTime;
        this.endTime = endTime;
        this.branchMode = branchMode;
        this.branchIndex = -1; // default -1 means no index assigned
        this.leftIndex = -1; // default -1 means no children
        this.rightIndex = -1;
        this.nodeType = -1; // default -1 means no type assigned
    }
}
