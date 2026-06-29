import numpy as np

import pandas as pd


REGIONS = ["SHA", "SGP", "FRA", "IAD", "NRT", "SYD"]
EVENT_TYPES = [
    "purchase",
    "refund",
    "click",
    "view",
    "search",
    "add_cart",
    "checkout",
    "share",
]


class DataGenerate:
    params = [
        [1_000_000, 5_000_000],
        [1_000, 50_000],
        [100, 200],
    ]
    param_names = ["rows", "users", "categories"]

    def setup(self, rows, users, categories):
        self.rng = np.random.default_rng(seed=42)
        self.user_id = self.rng.integers(0, users, size=rows)
        self.timestamp = pd.date_range("2025-01-01", periods=rows, freq="s")
        self.amount = self.rng.uniform(1.0, 500.0, size=rows)
        self.region = self.rng.choice(REGIONS, size=rows).astype(object)
        self.event_type = self.rng.choice(EVENT_TYPES, size=rows).astype(object)
        self.category_id = self.rng.integers(0, categories, size=rows)
        self.rating = self.rng.uniform(1.0, 5.0, size=rows)
        self.mask = self.rng.random(rows) < 0.15

        self.df = pd.DataFrame(
            {
                "user_id": self.user_id,
                "timestamp": self.timestamp,
                "amount": self.amount,
                "region": self.region,
                "event_type": self.event_type,
                "category_id": self.category_id,
                "rating": self.rating,
            }
        )

        self.df.loc[self.mask, "rating"] = np.nan


class Constructor:
    params = [
        [1_000_000, 5_000_000],
        [1_000, 50_000],
        [100, 200],
    ]
    param_names = ["rows", "users", "categories"]

    def setup(self, rows, users, categories):
        self.rng = np.random.default_rng(seed=42)
        self.user_id = self.rng.integers(0, users, size=rows)
        self.timestamp = pd.date_range("2025-01-01", periods=rows, freq="s")
        self.amount = self.rng.uniform(1.0, 500.0, size=rows)
        self.region = self.rng.choice(REGIONS, size=rows).astype(object)
        self.event_type = self.rng.choice(EVENT_TYPES, size=rows).astype(object)
        self.category_id = self.rng.integers(0, categories, size=rows)
        self.rating = self.rng.uniform(1.0, 5.0, size=rows)
        self.mask = self.rng.random(rows) < 0.15

    def time_construction(self, rows, users, categories):
        df = pd.DataFrame(
            {
                "user_id": self.user_id,
                "timestamp": self.timestamp,
                "amount": self.amount,
                "region": self.region,
                "event_type": self.event_type,
                "category_id": self.category_id,
                "rating": self.rating,
            }
        )

        df.loc[self.mask, "rating"] = np.nan


class GroupByAgg(DataGenerate):
    def time_groupby_only(self, rows, users, categories):
        self.df.groupby("user_id").ngroups

    def time_groupby_agg_mux(self, rows, users, categories):
        self.df.groupby("user_id").agg(
            total_amount=("amount", "sum"),
            avg_amount=("amount", "mean"),
            txn_count=("amount", "count"),
            unique_categories=("category_id", "nunique"),
            avg_rating=("rating", "mean"),
            region_count=("region", "nunique"),
        )

    def time_groupby_agg_sum(self, rows, users, categories):
        self.df.groupby("user_id").agg(
            total_amount=("amount", "sum"),
        )

    def time_groupby_agg_mean(self, rows, users, categories):
        self.df.groupby("user_id").agg(
            avg_amount=("amount", "mean"),
            avg_rating=("rating", "mean"),
        )

    def time_groupby_agg_count(self, rows, users, categories):
        self.df.groupby("user_id").agg(
            txn_count=("amount", "count"),
        )

    def time_groupby_agg_nunique(self, rows, users, categories):
        self.df.groupby("user_id").agg(
            unique_categories=("category_id", "nunique"),
            region_count=("region", "nunique"),
        )


class MergeJoin(DataGenerate):
    def setup(self, rows, users, categories):
        super().setup(rows, users, categories)

        self.user_dim = (
            self.df.groupby("user_id")
            .agg(
                first_event=("timestamp", "min"),
                last_event=("timestamp", "max"),
            )
            .reset_index()
        )
        self.user_dim["tenure_days"] = (
            self.user_dim["last_event"] - self.user_dim["first_event"]
        ).dt.total_seconds() / 86400.0

        cats = self.df["category_id"].unique()
        self.cat_dim = pd.DataFrame(
            {
                "category_id": cats,
                "category_name": [f"cat_{c}" for c in cats],
            }
        )

    def time_merge_user(self, rows, users, categories):
        self.df.merge(
            self.user_dim[["user_id", "tenure_days"]],
            on="user_id",
            how="left",
        )

    def time_merge_category(self, rows, users, categories):
        self.df.merge(
            self.cat_dim,
            on="category_id",
            how="left",
        )


class TimeSeries(DataGenerate):
    def setup(self, rows, users, categories):
        super().setup(rows, users, categories)
        self.hourly = self.df.set_index("timestamp").resample("h")["amount"].sum()

    def time_resample(self, rows, users, categories):
        self.df.set_index("timestamp").resample("h")["amount"].sum()

    def time_rolling_24h_mean(self, rows, users, categories):
        self.hourly.rolling(window=24, min_periods=1).mean()

    def time_rolling_24h_std(self, rows, users, categories):
        self.hourly.rolling(window=24, min_periods=1).std()


class StringCategorical(DataGenerate):
    def setup(self, rows, users, categories):
        super().setup(rows, users, categories)

        self.df_cat = self.df.copy()
        self.df_cat["region"] = self.df_cat["region"].astype("category")
        self.df_cat["event_type"] = self.df_cat["event_type"].astype("category")

    def time_to_categorical(self, rows, users, categories):
        df = self.df
        df["region_category"] = df["region"].astype("category")
        df["event_type_category"] = df["event_type"].astype("category")

    def time_string_ops(self, rows, users, categories):
        self.df["region"].str.upper()
        self.df["region"].str.contains("A", regex=False)
        self.df["event_type"].str.len()

    def time_value_counts(self, rows, users, categories):
        self.df_cat["region"].value_counts()
        self.df_cat["event_type"].value_counts()


class PivotCrosstab(DataGenerate):
    def time_pivot_table(self, rows, users, categories):
        self.df.pivot_table(
            values="amount",
            index="region",
            columns="event_type",
            aggfunc="mean",
        )

    def time_crosstab(self, rows, users, categories):
        pd.crosstab(self.df["region"], self.df["event_type"])


class ApplyVectorised(DataGenerate):
    def setup(self, rows, users, categories):
        super().setup(rows, users, categories)

        self.df_sample = self.df.head(min(rows, 500_000)).copy()

    def time_apply_rowwise(self, rows, users, categories):
        def _score_row(row):
            base = row["amount"] * 0.8
            if row["event_type"] == "purchase":
                base *= 1.5
            elif row["event_type"] == "refund":
                base *= 0.3
            if pd.notna(row["rating"]):
                base += row["rating"] * 10
            return base

        df = self.df_sample.copy()
        df.apply(_score_row, axis=1)

    def time_vectorised(self, rows, users, categories):
        df = self.df_sample

        multiplier = pd.Series(1.0, index=df.index)
        multiplier = multiplier.where(df["event_type"] != "purchase", 1.5)
        multiplier = multiplier.where(df["event_type"] != "refund", 0.3)
        rating_bonus = df["rating"].fillna(0.0) * 10
        df["score_vec"] = df["amount"] * 0.8 * multiplier + rating_bonus
