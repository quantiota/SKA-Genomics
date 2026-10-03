"""Test QuestDB connection and table operations."""
import psycopg2
from config import QDB_CONFIG


def test_questdb_connection():
    """Test basic QuestDB connection and operations."""
    try:
        conn = psycopg2.connect(**QDB_CONFIG)

        with conn.cursor() as cur:
            cur.execute("SELECT 1;")
            print(f"✅ QuestDB connection successful: {cur.fetchone()}")

            cur.execute("SHOW TABLES;")
            print(f"📋 Existing tables: {[t[0] for t in cur.fetchall()]}")

            cur.execute("DROP TABLE IF EXISTS test_genome;")
            cur.execute("""
                CREATE TABLE test_genome (
                    record_id SYMBOL,
                    value DOUBLE,
                    timestamp TIMESTAMP
                ) TIMESTAMP(timestamp);
            """)
            print("✅ Test table creation successful")

            cur.execute("DROP TABLE test_genome;")
            print("✅ Test table cleanup successful")

        conn.close()
        print("🎯 QuestDB validation PASSED")

    except Exception as e:
        print(f"❌ QuestDB connection failed: {e}")
        return False
    return True


if __name__ == "__main__":
    test_questdb_connection()
