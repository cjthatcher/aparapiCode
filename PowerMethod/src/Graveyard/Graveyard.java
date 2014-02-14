package Graveyard;

import java.sql.Connection;
import java.sql.SQLException;
import java.sql.Statement;

public class Graveyard {

	private static void logData(Connection connection, long scenario_id, int method, long duration, int size, int range, int result) throws SQLException
	{
		String sql = "INSERT INTO knapsack_data (method, scenario_id, runtime, result, size, range) VALUES (" + method + ", " + scenario_id + ", " + duration + ", " + result + ", " + size + ", " + range + ")";
		Statement s = connection.createStatement();
		s.executeUpdate(sql);
	}
}
