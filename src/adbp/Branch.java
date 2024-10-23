package adbp;

// class to store all information about a branch
public class Branch {
    public int startNode;
    public int endNode;
    public double startTime;
    public double endTime;
    public String branchType;
    public int branchIndex;
    public int leftIndex;
    public int rightIndex;

    public Branch(int startNode, int endNode, double startTime, double endTime, String branchType) {
        this.startNode = startNode;
        this.endNode = endNode;
        this.startTime = startTime;
        this.endTime = endTime;
        this.branchType = branchType;
        this.branchIndex = -1; // // default -1 means no index assigned
        this.leftIndex = -1; // default -1 means no children
        this.rightIndex = -1;
    }

    public void setIndex(int branchIndex) {
        this.branchIndex = branchIndex;
    }

    public void setLeftIndex(int leftIndex) {
        this.leftIndex = leftIndex;
    }

    public void setRightIndex(int rightIndex) {
        this.rightIndex = rightIndex;
    }
}
