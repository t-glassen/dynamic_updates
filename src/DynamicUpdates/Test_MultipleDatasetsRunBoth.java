package projects.DynamicUpdates;

import java.io.IOException;

public class Test_MultipleDatasetsRunBoth {

	public static void main(String[] args) throws IOException {
		Test_MultipleDatasetsWithIncreasingNCCPivot.main(null);
		Test_MultipleDatasetsWithMultipleAlphaCCPivot.main(null);
	}

}
